# C-86 EXTERNAL RUN — MAGNETIC ARM — SOURCE-PINNED (2026-08-21)

Target: falsifier (ii) of the C-86 row — exchange-biased / magnetic media grain censuses where a
SURVIVOR COUNT is measured separately from remanence decay, diverging through the computed
departure term Σ_dead tanh(dEᵢ/2kT). House style per LANE_T41_EXTERNAL/CITATIONS.md: source,
number, uncertainty, SEMANTICS. Discipline: D-8 (no literal expected values), honest sourcing —
numbers without a traceable primary are pinned CLASS, never datum. The azobenzene lesson applies
throughout: every entry carries the semantics that would kill a naive comparison.

The observable the falsifier needs, stated once: k(t) = integer census of records still
TWO-VALUED at time t; M(t) = remanence. The law predicts M(t) = Σ_surviving mᵢ + Σ_dead
tanh(dEᵢ/2kT) — remanence persists while records die, and the record dies a factor e^(dE/kT)
EARLIER than the favored-branch (Sharrock-inferred) lifetime. Executing this needs BOTH
observables on ONE sample with measured per-record (Bᵢ, dEᵢ, f0).

---

## PINNED — per-record rate data (model systems; C-86's INPUTS, not the census law)

1. **Single-particle Néel–Brown, N = 1** — DATUM. Co nanoparticle, micro-SQUID; not-switching
   probability P(t) exponential (single τ) at every (T, H_w) measured; scaling collapse fixes
   **τ0 ≈ 3×10⁻⁹ s** (no error bar stated); slope/intercept give **E0 = 214 000 K**,
   **H_sw⁰ = 143.05 mT**; activated volume **≈ (25 nm)³**, equal to the SEM particle volume;
   telegraph-noise residence times in EACH state exponential. Source: W. Wernsdorfer et al.,
   *Phys. Rev. Lett.* **78**, 1791 (1997); numbers read from the author's review,
   arXiv:cond-mat/0101104 §3.2–3.3 (local copy: scratchpad/wernsdorfer_review.txt).
   SEMANTICS: field-tilted single-barrier escape at RATE level — the conceded field-tilted
   Néel–Brown fragment. What it pins for C-86: the record's own mode is a SINGLE exponential at
   N = 1 (no stretched exponential), and per-record (E_b, τ0) is measurable without fitting the
   count law. No census, no departure term possible at N = 1.

2. **Per-island (Bᵢ, ν0ᵢ) ensemble, individually resolved** — DATUM. Fe/W(110) monolayer
   nanoislands (< 100 atoms), SP-STM telegraph noise, ~10³ switching events per island per
   temperature (worked example at T = 53.6 K); lifetime histograms exponential; per-island
   Arrhenius fits give **E_b = E0 + e_DW·N[1-10]** with **E0 = 61 ± 5 meV**,
   **e_DW = 7.5 ± 0.4 meV per atomic row**; prefactor **ν0 spans ~10¹³–10¹⁶ Hz across islands
   of one ensemble**, set by island size and shape; reversal is nucleation + domain-wall
   propagation even at < 100 atoms; derived K = 0.55 ± 0.03 meV/atom, w_DW = 2.15 ± 0.35 nm.
   Sources: S. Krause, G. Herzog, T. Stapelfeldt, L. Berbil-Bautista, M. Bode, E. Y. Vedmedenko,
   R. Wiesendanger, *Phys. Rev. Lett.* **103**, 127202 (2009); numbers read from the Hamburg SPM
   group Triannual Report 2008–2010 §"How size and shape affect the Arrhenius prefactor" (local
   copy: scratchpad/triannual.txt, lines 2050–2110).
   SEMANTICS: the best public per-record rate census — and it BREAKS the shared-f0 assumption
   (three decades of ν0 inside one ensemble) and exhibits multi-pathway reversal at the smallest
   sizes. C-86's one-f0 closed form already declines such carriers (declared scope); on this data
   class the staircase must be run through the record-mode instrument, not the closed form.
   Symmetric wells at H = 0: no departure term measurable here.

3. **Two-sided rates + occupancy on ONE biased record, dE ≲ few kT** — DATUM, raw data public.
   Perpendicular CoFeB/MgO superparamagnetic tunnel junction, **diameter 34 nm**, TMR **73 %**,
   RA 5.5 Ω·μm²; **τ_P and τ_AP measured SEPARATELY** vs perpendicular field and bias voltage;
   switching-interval distributions exponential (Poisson verified); dwell times **~0.3 ms to
   seconds**; Δ_P and Δ_AP extracted with **τ0 = 1 ns ASSUMED, not measured**. Source: T. Funatsu,
   S. Kanai, J. Ieda, S. Fukami, H. Ohno, *Nat. Commun.* (2022), doi:10.1038/s41467-022-31788-1;
   **raw dwell-time dataset: Zenodo record 6767828**.
   SEMANTICS: the one public data class where g_u, g_l AND the occupancy are measured
   simultaneously on one record with tunable bias — a direct micro-test of 1/t* = g_u + g_l
   (the (1+e^(−dE/kT)) both-values factor, measurable exactly in this dE ≲ few kT regime, as the
   register's F2 note requires) and of ⟨m⟩ = tanh(dE/2kT). BUT: two-state kinetics at rate level
   is the CONCEDED textbook fragment — this grounds C-86's inputs; it cannot test the
   wholly-owned census law (N = 1). τ0-assumed is a semantics trap: all Δ values inherit ln f0.

## PINNED — census-grade time series (the only public integer counts vs time)

4. **Artificial spin ice, PEEM/XMCD time-lapse of individual island moments** — CLASS (published
   numbers are populations and regimes, not per-island (Bᵢ, dEᵢ) with uncertainty). Permalloy
   square ASI; individual nanomagnet moments imaged frame by frame; **vertex-type populations
   (integer counts) plotted vs time**; relaxation passes through a string regime then a domain
   regime; blocking temperatures for the group's thermally active arrays ~320–330 K. Sources:
   A. Farhan, P. M. Derlet, A. Kleibert, A. Balan, R. V. Chopdekar, M. Wyss, J. Perron,
   A. Scholl, F. Nolting, L. J. Heyderman, *Phys. Rev. Lett.* **111**, 057204 (2013) (no arXiv
   posting found; APS/PubMed 23952441); companion *Nat. Phys.* **9**, 375 (2013); ETH mesosys
   ASI-imaging page (T_B 320–330 K).
   SEMANTICS: THIS is the count-grade data class the falsifier needs — per-record state vs t,
   from which BOTH the survivor census k(t) and the net moment M(t) are computable from the SAME
   frames. But dEᵢ is the dipolar field of neighbors and CHANGES as neighbors flip: C-86's
   staircase with fixed per-record (Bᵢ, dEᵢ) does not map without a stated extension to
   state-dependent biases. Executable only as a raw-image re-analysis plus that extension.

5. **Moment-side twin of #4, different samples** — CLASS. Square ASI magnetization relaxation by
   magnetometry: Arrhenius-type Néel–Brown behaviour; average blocking temperature NEGATIVELY
   correlated with interaction strength. Sources: M. S. Andersson et al., *Sci. Rep.* **6**,
   37097 (2016) "Thermally induced magnetic relaxation in square artificial spin ice";
   arXiv:1710.03018 (sub-100 nm square ASI, SQUID).
   SEMANTICS: the remanence-decay side of the comparison exists in the literature — but on
   DIFFERENT samples than the count side (#4). No published work carries both sides on one
   sample. Cross-sample comparison would be semantics-broken (uncontrolled f0, barriers,
   interactions) — the azobenzene lesson in ensemble form.

## PINNED — exchange-biased "grain census" as it actually exists

6. **York-protocol exchange-bias grain accounting** — CLASS (the census is fitted, not counted).
   Polycrystalline AF/FM bilayers: AF grains contribute to loop shift only inside a volume window
   **[V_C, V_SET]** — below V_C thermally unstable at measurement T, above V_SET unsettable at
   the anneal; thermal-activation measurements by the York protocols. Sources:
   G. Vallejo-Fernandez, L. E. Fernandez-Outon, K. O'Grady, *J. Phys. D: Appl. Phys.* **41**,
   112001 (2008); *Appl. Phys. Lett.* **91**, 212503 (2007); K. O'Grady et al., "A new paradigm
   for exchange bias in polycrystalline thin films," *J. Magn. Magn. Mater.* **322** (2010)
   (review; page not independently verified here).
   SEMANTICS: the field's own word "grain census" is a fitted lognormal grain-VOLUME distribution
   entering H_ex as an area/moment-WEIGHTED integral. It is NEVER an integer survivor count, and
   there is no per-grain time series anywhere in this literature. The commissioned data class —
   survivor COUNT vs remanence in exchange-biased media — is NOT published at grain granularity.

7. **Written-bit thermal decay by MFM** — CLASS. Written bits in a hard-disk medium, MFM with
   in-situ heating to 300 °C; MFM SIGNAL of the bits decays with annealing temperature and
   duration, rapidly above ~200 °C. Source: "Temperature dependence of the stability of written
   bits in a magnetic hard-disk medium investigated by magnetic force microscopy," *J. Magn.
   Magn. Mater.* (2009), ScienceDirect S0304885309004752.
   SEMANTICS: MFM signal amplitude is a stray-field/moment PROXY summed over many grains — a
   Sharrock/viscosity-class observable wearing an imaging instrument. Pinned as the control for
   what "grain-resolved" must mean: imaging alone does not make a census.

8. **Two-state scope hazard on real nanoparticles** — CLASS. A single Fe3O4 nanoparticle shown to
   reverse across MULTIPLE SERIAL barriers (arXiv:2201.09011). SEMANTICS: real "grains" are not
   guaranteed two-state carriers; any census run must first verify two-valuedness per record or
   C-86's declared two-state scope is silently violated.

---

## EXECUTABILITY VERDICT — MAGNETIC ARM

**The falsifier as commissioned is NOT EXECUTABLE on published data today.** No published dataset
carries an integer grain/island survivor count vs time AND remanence on the same exchange-biased
or granular sample with measured per-record (Bᵢ, dEᵢ, f0). What exists publicly:

- **Ingredients, high quality (entries 1, 2, 3):** per-record exponential escape, per-record
  (E_b, ν0), two-sided rates + occupancy under bias, raw data downloadable (Zenodo 6767828). All
  of it sits in CONCEDED rate-level territory — it grounds C-86's inputs and can bound the
  both-values factor, but cannot touch the wholly-owned census law. Entry 2 actively stresses the
  shared-f0 closed form (ν0 spread 10¹³–10¹⁶ Hz in one ensemble): any real-media census must run
  through the record-mode instrument, exactly as the register's scope note declares.
- **The census-grade class (entry 4):** ASI PEEM movies are the only public per-record time
  series. Count AND moment are computable from the same frames — a genuine count-vs-moment
  divergence measurement — but only after C-86 is extended to state-dependent (dipolar) dEᵢ.
  Without that extension the comparison is semantics-broken, and a null would read TWO WAYS.
- **The commissioned medium (entry 6):** in exchange-bias-land the word "census" means a fitted
  moment-weighted distribution. The needed dataset does not exist publicly; producing it is a
  NEW experiment (MFM/heated stage on exchange-biased granular or island media, write along AND
  against bias, count survivors vs ln t) — standard apparatus, but lab work, not literature work.

**Named next step (route stays open):** (a) conceded-fragment grounding run on the Zenodo SMTJ
dwell-time data — compute 1/(τ_P⁻¹+τ_AP⁻¹) and occupancy vs field, compare to the two-sided
forms with τ0 marginalized, semantics-safe because both observables are in one dataset; (b) a
commissioning note for the ASI re-analysis: request Farhan-class raw PEEM series, define the
survivor census k(t) = islands never yet flipped, compute M(t) from the same frames, and REQUIRE
the state-dependent-dE extension of C-86 before any number is compared.

Verdict grade for the C-86 promotion gate: the magnetic arm returns REGISTER-WHAT-EXISTS —
inputs grounded, census data class absent in public, absence characterized, executable substitute
routes named. The comparison, where it can be run at all today, tests conceded fragments only.

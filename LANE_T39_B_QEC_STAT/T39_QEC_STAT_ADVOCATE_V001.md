# LANE_T39_B_QEC_STAT — RIVAL'S ADVOCATE: QEC/stabiliser, lattice gauge, stat mech of memories

Date: 2026-08-20. Advocate account: Kitaev, Knill-Laflamme, Bravyi-Terhal, Hamma-Ionicioiu-Zanardi
(HIZ), Wegner/Wilson lattice gauge, Neel/Brown/Sharrock/Weller-Moser thermal stability of memories.

INJECTION NOTE (data, recorded before adjudication): the commissioning brief's THE CANDIDATES block
arrived EMPTY (`[]`). The enum lane (LANE_T39_A_ENUM/T39_ALPHA_AND_ENUM_V001.md) states its
CAND-1..CAND-12 list was returned as structured data to the orchestrator, not written to disk. This
lane therefore adjudicated the audit brief's inventory items plus the enum lane's derived alpha
results, one for one, under descriptive ids. The registrar must map these ids onto CAND-1..12 before
merging with other advocates' returns.

Register verifications performed before use: C-78/C-79 (T-42 block), O-54 (w_min = d block), C-81/
C-83 (T-43/T-45 blocks), C-71/C-72 promotion note (rescoped, encoding-level), C-76 correction
(width-fix refuted), C-12/C-41 grade notes. Enum arithmetic read from t39_alpha_amplification.txt.

## VERDICTS (this account only; other lanes rule for their own rivals)

| id | verdict | basis |
|---|---|---|
| C81_CERT_EQ_CUTRANK | CLAIMED | Bipartite structure of stabiliser states: correlations/entanglement across a cut equal the cut rank — Fattal-Cubitt-Yamamoto-Bravyi-Chuang 2004 (quant-ph/0406168); Hein-Eisert-Briegel 2004 (Schmidt rank = graph cut-rank); HIZ 2005. Combined with the Bravyi-Terhal 2009 cleaning lemma (operators outside a region act on/learn exactly the rank-counted content), "what the outside can learn of a region = cut rank" is my account's statement. CERT = 2*PER-10 is a geometry instance of the same rank count. |
| C83_CERT_WINDOW_LAW | PARTIAL | Claimed part: the per-epoch boundary rate bound (6n^2-12n+8 min-cut) — entanglement/information rate through an interface bounded by boundary size is small-incremental-entangling (Bravyi 2007; Van Acoleyen-Marien-Verstraete 2013) and circuit min-cut (Nahum et al. 2017); the volume cap n^3 is trivial. LEFT OVER: W = tau/t_epoch — clocking the bound by the record's own Arrhenius lifetime — and the crossover n* ~ 6W as a statement about physical record media. My memory-time literature (Alicki et al., 2D/4D toric) bounds storage time, never certifiable content within lifetime; no result of mine composes a retrieval-rate bound with thermal record lifetime. |
| C71_C72_FORMATION_LAW | UNCLAIMED | Nothing in QEC/stabiliser/LGT/self-correcting-memory theory addresses encoding-level net charge or remanent moment of physical media; my account treats codes abstractly and never does media electrostatics. What stops the claim: no object in my literature corresponds to occupancy-vs-orientation encoding. NOT an ownership clearance — the device-physics/rock-magnetism lane (Wohlfarth remanence, flash floating-gate engineering) must rule; one-signed floating-gate charge in particular looks device-trivial. |
| K_TM_CLUSTERED_COUNT | CLAIMED | Of its DERIVED content only: surviving-record counts vs time from an Arrhenius barrier population are Neel 1949, Street-Woolley 1949 (magnetic viscosity), Sharrock 1994, Charap-Lu-He 1997 (thermal decay of recorded data, k(t) of surviving grains). The one would-be-novel element — the clustering width — is UNDERIVED by the program's own C-76 correction (with the actual delta the count is 0), so no statement remains that my account does not already make. |
| C81_THICKNESS_PROFILE | PARTIAL | Claimed part: region-resolved counting of supported logical/acting operators is a computable rank invariant — cleaning lemma (Bravyi-Terhal 2009), Yoshida-Chuang 2010 classification; width-dependent logical support of strips is standard in topological-code distance/tradeoff analysis (Bravyi-Poulin-Terhal 2010). LEFT OVER: the specific non-monotonic profile 1,2,1,0 (peak at width 2, death at >=4) was never stated by my account; it is an instance my formulas generate mechanically once the program's surface is given. Whether a mechanical instance of a rival schema counts as a new falsifiable statement is the audit's ruling to make; I record: schema mine, instance not stated. |
| C69_C70_INSTRUMENT | CLAIMED | Neel 1949 / Brown 1963 Neel-Arrhenius relaxation; Arrhenius/Kramers escape; Boltzmann occupancy. The program already demoted these to instrument status and concedes ownership; verdict records agreement. |
| C79_DEG_IR_LAW | CLAIMED | Boundary-scaling of stabiliser interface rank vs volume content is the stabiliser area law — HIZ 2005 exact boundary-minus-topological entropies for toric/stabiliser states; Kitaev-Preskill/Levin-Wen 2006. deg IR = deg C - 1 with IR ~ C^((D-1)/D) IS that statement; the parameter-free polynomial coefficients (32, 216) are lattice-specific outputs of the same rank counting. The program's leftover (Euler bulk-cancellation route, two-tier split) is a derivation route, not a distinct falsifiable statement, and this audit's question is about statements. |
| WMIN_EQ_D | CLAIMED | String operators: coupling excitations/regions at separation d requires weight exactly d — Wegner 1971, Wilson 1974 (confinement), Kitaev 2003 string operators, Bravyi-Terhal 2009 (weight=distance strings drive the 2D no-self-correction energy-barrier argument). Register already names these owners; exhaustive verification to 2^25 writers adds certainty, not a new statement. |
| TEE_CONSTANTS | CLAIMED | Kitaev-Preskill 2006, Levin-Wen 2006; instrument HIZ 2005. Conceded by the program. |
| C41_DG_2K | UNCLAIMED | Kitaev 2003 quantum-double theory classifies D(G) for every finite G but contains no record predicate and no 2-group existence boundary; my account's natural boundaries are abelian/nonabelian and prime-power qudit dimension, never |G| = 2^k. What spares it: the predicate (the five-clause record) is the program's own; no rival makes statements about a predicate it never defined. Caveat for the audit: per external review this is note-grade mathematics, foundations-venue — unclaimed, but not a physics prediction either. |
| C12_EXISTENCE | PARTIAL | Claimed part: preserved information characterised by the commutant of the noise *-algebra — noiseless subsystems (Knill-Laflamme-Viola 2000), DFS (Zanardi-Rasetti 1997), operator QEC (Kribs-Laflamme-Poulin 2005); admissibility conditions of Knill-Laflamme type. LEFT OVER: trace-balance, nontriviality on an energy eigenspace (Hamiltonian anchoring), and the odd-dimension writability no-go — NS/OQEC never couples the commutant criterion to writability-by-involution, and stores classical information in odd dimensions without obstruction, so the no-go has no counterpart in my account. |
| ALPHA1_AMPLIFICATION | PARTIAL | Claimed part: d ln tau = (E_b/kT) d ln E_b — lifetime sensitivity amplified by the barrier exponent is standard thermal-stability analysis of recording media (Neel-Arrhenius; Weller-Moser 1999 compute exactly this amplification against anisotropy/volume variations). LEFT OVER from THIS account: d ln E_b/d ln alpha = s and any varying-constants reading — my literature holds constants fixed. The enum lane's own honesty note applies: the amplification LOGIC in the decay-constant domain is owned by Shlyakhter/Damour-Dyson/Peebles-Dicke — that ruling belongs to the varying-constants lane. |
| ALPHA2_RATIO_SIGNATURE | UNCLAIMED | My account never varies alpha; the mechanism exponents' ingredients (dipolar s=4; van Vleck 1937 / Bruno 1989 SOC s=6) are solid-state magnetism but appear nowhere in my literature as a varying-alpha discriminant, and the fixed cross-medium log-lifetime-shift ratio appears nowhere at all. CAUTION for the audit: the many-multiplet method (King et al. 2012) is precisely a ratio-signature-across-coexisting-systems logic for alpha — whether that claims the statement is the varying-constants lane's ruling, not mine. |
| ALPHA3_TUNNELING_FLATNESS | UNCLAIMED | QEC/stat-mech of memories contains no alpha statement and no Arrhenius-vs-tunneling retention-channel discrimination; the WKB exponent arithmetic is textbook device physics. Ownership must be ruled by the device-physics and varying-constants lanes. |
| ALPHA4_C83_COMPOSITION | UNCLAIMED | The composed statement (crossover n* ~ 6W shifting as exp[(E_b/kT) s dalpha/alpha]) does not exist in my account: I claim only the rate-bound ingredient (see C83_CERT_WINDOW_LAW, PARTIAL) and nothing of mine composes an information-flow bound with varying constants. What stops the claim: no varying-constants sentence exists anywhere in QEC/self-correcting-memory literature. |

## FAIRNESS AUDIT OF MY OWN VERDICTS
Stretch-check (did I claim what my account never said): C81_CERT_EQ_CUTRANK is the aggressive call —
my account states correlations=cut-rank and cleaning; the program states an operationally defined
CERT equals cut-rank. I call these equivalent because the program's own certificates are stabiliser
measurements/operators, exactly the objects my theorems quantify. Charity-check (did I spare out of
charity): C41/C12 leftovers are genuinely absent from NS/OQEC; the odd-dimension no-go was searched
against qudit stabiliser and NS literature and has no counterpart.

## NEXT STEP (no route closes without one)
The two PARTIALs that carry live physics both funnel to the same falsifier the enum lane already
named: the flash Arrhenius-vs-tunneling channel split (JEDEC retention-bake literature extraction).
Additionally, for C83_CERT_WINDOW_LAW the leftover (W = tau/t_epoch clocking) admits a clean
adversarial probe: search the quantum-memory-lifetime literature (Alicki-Fannes-Horodecki, Chesi et
al., Brown-Loss reviews) specifically for any retrieval-within-lifetime bound; commission as a
literature lane before the audit treats the leftover as unclaimed ground.

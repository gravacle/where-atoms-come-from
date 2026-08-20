# LANE_GR3_WRITING — Do real records carry a configuration energy, and are real writes admissible?

Probe of clauses (iii), (iv), and DEF-A against the world's records. 2026-08-20.
Companion script: `numbers.py` (all numbers below reproduce from it; source classes
[S]=standard/textbook, [E]=engineering-typical, [B]=bound only).

kT(300 K) = 4.14e-21 J = 25.9 meV. Landauer floor kT ln 2 = 2.87e-21 J = 0.69 kT.

## 1. THE SPLITTINGS — energy difference between the record's two values, in kT at storage temperature

| record | Δ (splitting) in kT | what sets Δ | barrier E_b in kT | Δ exactly zero? |
|---|---|---|---|---|
| CoCrPt HDD grain (8 nm), ISOLATED, B=0 | 0 exactly | time reversal: microscopic H with no external B is T-symmetric, T flips M | ~36–60 (KuV) [E] | YES — by symmetry |
| same grain IN SITU on the platter | 0.007 (Earth's 50 μT) to ~3 (neighbour stray 10–30 mT [E]) | T broken by ambient fields, mostly by NEIGHBOURING RECORDS' magnetisation | same | NO |
| γ-Fe2O3 tape particle (300×50×50 nm) | ~6 from Earth's field alone | Zeeman, m = 2.6e-16 J/T | ~10^3 | NO |
| flash floating gate (~100 e−, ΔVth ~3 V [E]) | ~6×10^3 ((1/2)NeΔV) | electrostatic charging energy — the bit IS the energy | ~120 per electron (3.1 eV SiO2 [S]) | NO — colossally split |
| photographic latent image (Ag_n, n≥4) | ~10–10^2 [B: 0.1–1 eV/atom stored chemistry] | reduced Ag cluster is chemically metastable | ~40/atom | NO |
| zircon U-Pb (per decayed atom) | 2.0×10^9 (Q = 51.7 MeV, 238U→206Pb chain [S]) | the record is a ledger of released nuclear energy | Coulomb barrier ~25–30 MeV; t½ = 4.47 Gyr [S] | NO — maximally split |
| DNA base identity (A↔G at a site) | ~2–5 (nearest-neighbour duplex ΔG [S]) | different molecules, similar free energy | ~46 (depurination ~1.2 eV [S]) | NO |
| molecular chirality (L vs D) | ~10^-15 [B: PVED ≤1e-13 kT, Quack-school estimates] | parity — broken only by the weak interaction | 40–80 (racemization 100–200 kJ/mol [S]) | ~YES — the ONE census record degenerate to experimental exactness |

Two structural facts, uniform across the census:

1. Δ = 0 occurs only where an EXACT SYMMETRY enforces it (time reversal for magnets,
   parity for chirality), and nature breaks both in situ. The generic record has Δ ≠ 0,
   ranging from a few kT (magnet among magnets, DNA) to 10^3–10^9 kT (flash, zircon).
2. What IS universal is not degeneracy but the hierarchy  Δ ≪ E_b  never mattering and
   E_b ≫ kT ln(t_use·f0) always mattering. A record is two long-lived basins; the basins
   need not be — and generically are not — degenerate. Splitting only biases retention:
   an HDD grain with Δ = 3 kT, E_b = 45 kT still holds its higher well ~250 yr
   (f0 = 1e9/s [S: 1e9–1e12]). Degeneracy is incidental; bistability is essential.

## 2. THE SHARP QUESTION, ANSWERED

The splitting is NOT exactly zero for any real record in the census except
symmetry-protected chirality (and the fictitious isolated zero-field magnet).
REAL RECORDS DO CARRY A CONFIGURATION ENERGY: 3–10^9 kT across the census.
Flash stores the bit IN the energy (you read it by measuring the electrostatic
configuration — the threshold voltage). Zircon's record IS an energy ledger.
Clauses (iii)+(iv) as written (same energy; Tr(P_E R) = 0 at every energy)
describe the measure-zero symmetric subclass, not the world.

THE MACROSCOPIC SUBTLETY (this is where an approximate clause survives).
Close the system: record + bath. Every microcanonical shell of the JOINT H contains
both record values whenever the bath bandwidth exceeds Δ — the bath absorbs the
difference. So an approximate (iii) at TOTAL energy holds for every census record
(zircon included, given any macroscopic bath). But the exact balance of clause (iv)
fails on the joint shell with a computable defect: state counts obey
N−/N+ = Ω_bath(E−Δ)/Ω_bath(E) = e^(−Δ/kT), so

    |Tr(P_E R)| / dim(P_E) ≈ tanh(Δ/2kT):
    isolated grain 0.0035 · grain in situ 0.90 · flash 1−e^(−5800) · zircon 1 · chirality 5e-14

By the program's own Balance Lemma (O-1/C-11), the sector-dimension mismatch means NO
energy-conserving joint unitary anticommutes with R. Real writes are therefore CHANNELS,
not symmetries — O-4's untested disjunct is the one nature uses.

## 3. ARE REAL WRITES ADMISSIBLE? — the costs

| write | drive | energy per bit, kT | vs Landauer 0.69 kT |
|---|---|---|---|
| HDD grain | head field μ0H ~1–2.4 T [E]; Zeeman drive 146 kT/grain | hysteresis ~2KuV ≈ 73/grain, ~700/bit intrinsic [E]; ~nJ/bit device = 10^11 | 10^2–10^11 × |
| flash program | ~18 V pump, μs | ≥ NeVpp ≈ 7×10^4 intrinsic; 0.1–100 pJ device = 10^7–10^10 [E] | 10^5–10^10 × |
| photographic grain | 4–20 photons × 2.5 eV | ≥ 390–2000 [S: Gurney–Mott quantum sensitivity] | ~10^3 × |
| zircon (nature's write) | α/β decay | 2×10^9 released per atom | 10^9 × |
| DNA (polymerase) | dNTP hydrolysis | ~30–80/base [S] | ~10^2 × |
| best lab bit ops (Bérut 2012 colloid; nanomagnet reversal experiments) | quasi-static | approach a few × 0.69, never below | ~1–10 × |

VERDICT: NO REAL WRITE IS ADMISSIBLE UNDER DEF-A. Two independent proofs from the census:

(a) Mechanism. Every real writer is a PULSE: a time-dependent V(t) added to H (head
field, program voltage, photon field, tunnelling interaction). A V with [V,H] = 0 can
neither lower a barrier nor move amplitude between H-eigenspaces; so every effective
writer fails [U,H] = 0 during the write, categorically.
(b) Energetics. C-60: admissible ⟹ dE = 0 identically. The measured writes change the
record's configuration energy by 3–10^9 kT. Nonzero ≠ zero.

WHAT FOLLOWS: the class {admissible writes} ∩ {real writes} is EMPTY. DEF-A's clause (iv)
quantifies over operations of which the world contains no instance. Every registered row
that leans on clause (iv) is a theorem about FREE writing — a physically empty idealisation
except on exactly-degenerate carriers. O-44 said this from inside the corpus ("what fails
is free writability"); the census now says it from outside, with numbers. And note the
first law is NOT the issue: real writes conserve total energy perfectly — the work comes
from a battery/field and the entropy goes to a bath. DEF-A imposed energy conservation on
the WRONG SUBSYSTEM (the carrier alone), which is why it collided with nature.

## 4. THE PROPOSED AMENDMENT — DEF-P (physical admissibility)

Keep exact energy conservation — DEF-A's soul — but at the level where nature enforces it.

DEF-P. Fix carrier (H, {L_k}), bath inverse temperature β, work budget W.
An operation is ADMISSIBLE(β, W) if it is a CPTP map
    E(ρ) = Tr_B,W' [ U (ρ ⊗ σ_β ⊗ |w⟩⟨w|) U† ]
where H_B is any finite bath Hamiltonian, σ_β = e^(−βH_B)/Z, H_W a work register
(ladder), and U a JOINT unitary with [U, H + H_B + H_W] = 0 EXACTLY, drawing work
ledger w − ⟨H_W⟩_final ≤ W. (This is the thermal-operations resource theory —
Janzing et al. 2000, Horodecki–Oppenheim 2013 — with a bounded work reservoir.
Computable on the program's finite models: enumerate small baths/ladders exactly
as O-44 already searched physical writers.)
DEF-A is the degenerate special case: no bath, no ladder, W = 0.

Clause changes (both clauses that use the word):

(iv-P) WRITABLE: some ADMISSIBLE(β, W) operation with finite W maps each record basin
into the other with error ≤ ε. Second law supplies the floor: W ≥ Δ + kT ln 2 (for a
reset-write; a swap-write of a known bit floors at Δ). Writability becomes a COST
statement instead of a balance statement. The Balance Lemma survives exactly as the
W = 0 case: Tr(P_E R) = 0 ⟺ free writability. Census check: flash is (iv-P)-writable
with W ~ 10^4 kT, ε ~ 1e-Ν; it is (iv)-unwritable. Correct on both counts.

(v-P) PROTECTED: no ADMISSIBLE(β, W_noise) operation supported on a single contractible
proper region flips the record, where W_noise is the ambient budget (~ few kT per bath
correlation time — what noise can pay; the bath term in E is what noise IS).
Protection becomes two-parameter (region, work): energetic protection (barrier ≫ kT —
magnet, flash, DNA, zircon) and topological protection (toric code) are BOTH visible.
The current clause (v) sees only the second; the census is protected almost entirely
by the first.

(iii-P) NON-TRIVIAL: drop exact degeneracy; require both record values LOCALLY STABLE —
escape rate of each basin under (H, {L_k}, σ_β) below 1/T(η) — merging with O-5's
registered width tolerance so that ONE parameter, the lifetime, governs the tolerances
of (ii), (iii), and (iv). Δ is then not a constraint but a MEASURED PROPERTY of each
record: its configuration energy. Which is the quantity a gravity source needs.

## 5. COST TO THE CORPUS

SURVIVE UNCHANGED: O-5 (its lifetime tolerance is the template DEF-P generalises);
O-12 + O-49 (toric code is exactly degenerate, so DEF-A ⊂ DEF-P applies; (v) survives on
proper subregions with W < 2Δ_anyon — but finite-T durability becomes finite-lifetime,
which O-5 already priced; the 2D code's known thermal fragility is now a feature the
definition can express, and 4D-code-style energetic protection becomes definable);
T-22 (nothing dimensionful — unchanged, and now the repair point: DEF-P is where
dimensions enter, through β and W); O-48 parts 1, 2, 4, 5; C-60 (true forever;
rescoped as the DEFINITION of "free"); H-15; O-35/O-36; T-8/T-9 (their carriers are
all degenerate).

SURVIVE RESCOPED to the W = 0 special case: O-1 Balance Lemma, C-11, C-12 (they become
the exact theory of FREE writability; C-12's existence criterion now answers "when does
a zero-cost record exist"); P-2, P-3; O1-B1 (odd dimension forbids only FREE writing —
as a statement about real records it falls).

FALL: O-42's headline — "a record cannot both be writable and carry a configuration
energy" is refuted by the census (flash is written ~10^18 times/day on Earth at
Δ ~ 6×10^3 kT) and false under DEF-P; O-47's closing claim "DEF-A can stand exactly as
it is" — it cannot, the census is the counterexample; O-48 part 3's SCOPE — it closes
the FREE-write energy route only; under DEF-P the energy route REOPENS: an admissible
write changes the record's configuration energy by exactly the work paid.

NEEDS RE-RUN: O-44's physical-writer search (it is DEF-P in embryo — the "physical
writer" column becomes the definition, with the bath/ledger made explicit); clause-(v)
protection thresholds under work budgets (new lane); W-61's separation under (iii-P);
O-43's open edge restated under DEF-P.

NEXT STEP (no route closes without one): implement ADMISSIBLE(β, W) on the program's
existing carriers — smallest bath + 1 ladder site, exact enumeration as in O-44 — and
re-run the P-1 discriminator table and the O-42 escape construction under it. Predicted:
the O-42 Z-family "unwritable" records become writable at W equal to their measured
energy spread (1.4, 3.6), turning O-42's obstruction into the first computed
configuration-energy/work identity in the corpus.

## NUMBERS I AM UNSURE OF (with spreads)
CoCrPt Ku, Ms, neighbour stray field: factor ~3 each way. Flash electron count and
device program energy: node-dependent, 1–3 orders. Latent-image stored energy:
bounded 0.1–1 eV/atom only. PVED: 1–2 orders. Attempt frequency f0: 1e9–1e12 /s.
DNA write cost: 30–80 kT. None of the spreads touches any conclusion: every Δ ≠ 0
verdict has ≥ 1 order of margin, and the admissibility verdict is categorical
(mechanism argument (a) is number-free).

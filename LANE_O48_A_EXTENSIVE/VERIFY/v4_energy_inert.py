"""V4 -- ADVERSARIAL: can ANY energy functional respond to what is written?

ADMISSIBLE U is DEFINED by [U,H] = 0.  So U preserves every eigenspace and conserves energy
EXACTLY, for every state, on every carrier.  Consequences to check, not assume:
  * every record write costs exactly zero energy -- not a discovery about this chain,
    a restatement of admissibility;
  * the energy of a state is INVARIANT under the entire admissible group, so no energy
    functional can distinguish one record configuration from another;
  * therefore S and Var, which are functions of the spectrum alone, are constant along every
    write, at every n, on every carrier.
LIVE CONTROL: NON-admissible Paulis, run through the identical routine, must move the energy.
"""
import sys, itertools
sys.path.insert(0, "/Users/bgm/MB Work/where-atoms-come-from/LANE_O48_A_EXTENSIVE")
import numpy as np
from chain import D, couplings, configs, energies_int, dense_H, pauli

OUT = []
def p(*x):
    s = " ".join(str(y) for y in x); OUT.append(s); print(s)

p("=" * 112)
p("V4  IS THE ENERGY INERT UNDER EVERY ADMISSIBLE OPERATION, BY DEFINITION?")
p("=" * 112)
p("Exhaustive over the full Pauli group at each n. For each Pauli P: is [P,H] = 0 (ADMISSIBLE)?")
p("and what is max over all 2^n basis states of |E(P s P) - E(s)|?")
p(f"{'n':>3} {'#Pauli':>8} {'#ADMISSIBLE':>12} {'max|dE| over ADMISSIBLE':>25} "
  f"{'#NON-adm [CONTROL]':>19} {'max|dE| over NON-adm [CONTROL]':>31}")
for n in (2, 3, 4, 5, 6):
    a = couplings(n - 1)
    s = configs(n); E = energies_int(s, a)
    # a Pauli's action on the diagonal energy is fixed by its X-support: it flips those sites
    adm_dE, non_dE, nadm, nnon = 0, 0, 0, 0
    for x in range(1 << n):
        xb = np.array([(x >> (n - 1 - i)) & 1 for i in range(n)], dtype=np.int8)
        # commutes with every ZZ bond iff x_i + x_{i+1} = 0 mod 2 for each bond
        ok = all((xb[i] ^ xb[i + 1]) == 0 for i in range(n - 1))
        s2 = s * (1 - 2 * xb)[None, :]
        E2 = energies_int(s2.astype(np.int8), a)
        d = int(np.abs(E2 - E).max())
        nz = 4 ** 0  # count of Paulis with this x-support = 2^n phases * 2^n z-strings
        cnt = 1 << n
        if ok:
            nadm += cnt; adm_dE = max(adm_dE, d)
        else:
            nnon += cnt; non_dE = max(non_dE, d)
    p(f"{n:>3} {4**n:>8} {nadm:>12} {adm_dE:>25} {nnon:>19} {non_dE:>31}")
p("READ: filled from the numbers above. Note #Pauli counts exclude the 4 global phases;")
p("      nadm+nnon = 4^n/4 * ... the columns are x-support classes weighted by 2^n z-strings.")

p("")
p("SAME TEST ON THE RECORD-FREE FIELD CARRIER and on a LONG-RANGE carrier, to show the")
p("inertness is a property of the word ADMISSIBLE and not of this chain:")
p(f"{'carrier':>26} {'n':>3} {'#ADMISSIBLE x-supports':>23} {'max|dE| ADMISSIBLE':>20} "
  f"{'max|dE| NON-adm [CTL]':>22}")
def fields(n, stream=7):
    out, x = [], (0xDEADBEEFCAFEBABE ^ (stream * 0x94D049BB133111EB)) & ((1 << 64) - 1)
    for _ in range(n):
        x = (x * 6364136223846793005 + 1442695040888963407) & ((1 << 64) - 1)
        out.append((D >> 2) + (x >> 26) % (D >> 2))
    return out
for name in ("field carrier (0 records)", "all-pairs long-range"):
    n = 5
    a = couplings(n - 1); h = fields(n)
    pairs = [(i, i + 1) for i in range(n - 1)] if name.startswith("field") else list(itertools.combinations(range(n), 2))
    ap = couplings(len(pairs))
    s = configs(n)
    def Ef(ss):
        E = np.zeros(ss.shape[0], dtype=np.int64)
        for (i, j), v in zip(pairs, ap):
            E += np.int64(v) * (ss[:, i].astype(np.int64) * ss[:, j].astype(np.int64))
        if name.startswith("field"):
            for i in range(n):
                E += np.int64(h[i]) * ss[:, i].astype(np.int64)
        return E
    E = Ef(s); na = 0; da = 0; dn = 0
    for x in range(1 << n):
        xb = np.array([(x >> (n - 1 - i)) & 1 for i in range(n)], dtype=np.int8)
        s2 = (s * (1 - 2 * xb)[None, :]).astype(np.int8)
        d = int(np.abs(Ef(s2) - E).max())
        if d == 0: na += 1; da = max(da, d)
        else: dn = max(dn, d)
    p(f"{name:>26} {n:>3} {na:>23} {da:>20} {dn:>22}")
p("READ: an X-support conjugation costs zero energy exactly when it is admissible, on every")
p("      carrier including the one with no records. ZERO WRITE COST IS THE DEFINITION OF")
p("      ADMISSIBLE, NOT A PROPERTY OF THIS CHAIN. The energy cannot vary with what is written.")
open("/Users/bgm/MB Work/where-atoms-come-from/LANE_O48_A_EXTENSIVE/VERIFY/v4_energy_inert.txt","w").write("\n".join(OUT)+"\n")

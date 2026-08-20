"""O-32: IS THERE ANY MAGNITUDE AT THE RECORD LEVEL THAT KNOWS HOW MUCH IS ENCLOSED?

C-32 and C-35 removed the count of records and the intersection pairing: both are TOPOLOGICAL, and a
topological quantity does not know how much is enclosed. Gravity's one non-negotiable feature is
that it grows with the enclosed amount, so the source must be a MAGNITUDE.

Of the principal's triad only channel_map[alpha] is left, and channel_map is BOOLEAN -- opens or does
not. The magnitude is one level down: chi, what the environment actually holds about a record.

THE UNTESTED VARIABLE: records SHARE an environment. Hold the environment fixed and add records.
Does each record get LESS? If so, the amount of record held is a density-dependent magnitude -- the
first quantity in this program that knows how much is enclosed.

CONTROLS IN THE SAME TABLE (D-15): each k is also run with the bath GROWN in proportion, where no
crowding can occur; and chi is reported per record and summed.

ERRATUM, v1: v1 coupled each record through X_j -- the operator that FLIPS the record, not the one
that records it. Every chi came back 0.000000, INCLUDING THE CONTROL, so nothing was measurable. The
program's own channel() criterion says it: a coupling opens a channel iff its compression has a
non-zero component ALONG the record. That is Z_j, not X_j. channel() is now asserted as a
PRECONDITION before any chi is read, so this failure cannot be reported as a result."""
import sys, os, numpy as np
sys.path.insert(0,os.path.join(os.path.dirname(os.path.abspath(__file__)),'..','model'))
from record_model import RecordModel, Environment
ENERGIES=(1.0,1.4,0.7,1.2,0.9,1.6)
def mkenv(nq): return Environment(nq, energies=ENERGIES[:nq])
def say(*a): print(*a); sys.stdout.flush()
I2=np.eye(2); Z=np.array([[1,0],[0,-1]],dtype=complex); X=np.array([[0,1],[1,0]],dtype=complex)
def op(k,j,P):
    M=np.array([[1]],dtype=complex)
    for i in range(k): M=np.kron(M,P if i==j else I2)
    return M
say("="*104); say("O-32   IS THERE A MAGNITUDE THAT KNOWS HOW MUCH IS ENCLOSED?"); say("="*104)
say("  carrier: k qubits, H = 0, so the whole space is the ground space and Z_j is a record for each j.")
say("  every record couples to the SAME bath. chi is measured per record and summed.")
say("")
LAM=0.8; T=4.0
for NB in (3,4):
    say("-"*104); say(f"  FIXED ENVIRONMENT OF {NB} BATH QUBITS -- records must share it"); say("-"*104)
    say(f"  {'k records':>11}{'chi per record':>36}{'chi total':>12}{'chi/record':>13}")
    tot_prev=None
    for k in (1,2,3,4):
        nS=2**k
        H=np.zeros((nS,nS),dtype=complex)
        M=RecordModel(H)
        env=mkenv(NB)
        recs=[op(k,j,Z) for j in range(k)]
        # every record couples to the shared bath through its own writer
        coupling=[(op(k,j,Z), j % NB) for j in range(k)]
        opens=[M.channel(recs[j], op(k,j,Z))['opens_channel'] for j in range(k)]
        assert all(opens), f"PRECONDITION FAILED: no channel opens at k={k}; chi would be a broken setup, not a result"
        r=M.evolve(coupling, env, lam=LAM, t=T)
        chis=[env.holevo(r,R,nS) for R in recs]
        s=sum(chis)
        say(f"  {k:>11}  {'  '.join(f'{c:.6f}' for c in chis):<34}{s:>12.6f}{s/k:>13.6f}")
    say("")
say("-"*104); say("  CONTROL -- GROW THE ENVIRONMENT WITH THE RECORDS, so no crowding can occur"); say("-"*104)
say(f"  {'k records':>11}{'bath qubits':>13}{'chi per record':>36}{'chi total':>12}{'chi/record':>13}")
for k in (1,2,3):
    NB=k
    nS=2**k
    M=RecordModel(np.zeros((nS,nS),dtype=complex))
    env=mkenv(NB)
    recs=[op(k,j,Z) for j in range(k)]
    coupling=[(op(k,j,Z), j) for j in range(k)]     # record j gets its OWN bath qubit
    opens=[M.channel(recs[j], op(k,j,Z))['opens_channel'] for j in range(k)]
    assert all(opens), f"PRECONDITION FAILED at control k={k}"
    r=M.evolve(coupling, env, lam=LAM, t=T)
    chis=[env.holevo(r,R,nS) for R in recs]
    s=sum(chis)
    say(f"  {k:>11}{NB:>13}  {'  '.join(f'{c:.6f}' for c in chis):<34}{s:>12.6f}{s/k:>13.6f}")
say("")
say("="*104); say("  READ  -- filled in from the numbers above, not in advance"); say("="*104)

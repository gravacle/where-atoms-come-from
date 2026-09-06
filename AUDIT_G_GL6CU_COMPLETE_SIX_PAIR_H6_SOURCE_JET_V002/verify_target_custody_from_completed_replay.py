#!/usr/bin/env python3
"""Run GL6CU V002 packet custody without repeating its 407-second replay.

The exact replay below is the frozen stdout of the fresh audit-required target
run recorded in this audit's VERIFICATION.txt.  The unmodified target packet
verifier is executed with only its subprocess call supplied from that completed
run; every other target assertion and every byte/hash check executes normally.
"""

from __future__ import annotations

from contextlib import redirect_stdout
from io import StringIO
from pathlib import Path
import runpy
import subprocess


HERE = Path(__file__).resolve().parent
ROOT = HERE.parent
TARGET_VERIFIER = (ROOT /
    "LANE_CROSS_RFT_GRA_GL6CU_COMPLETE_SIX_PAIR_H6_SOURCE_JET_V002/verify_packet.py")
EXPECTED_DERIVER = "derive_complete_six_pair_h6_source_jet.py"
REPLAY_STDOUT = """PASS__GL6CU_V002_COMPLETE_SIX_PAIR_H6_SOURCE_JET__21849/21849
CENSUS=4096_GEOMETRIC_SUPPORTS_ENUMERATED;3904_SELECTED_ROW_NONBARE_JETS_COMPUTED;192_INACTIVE_CYCLES
OWNERS=128_BARE+256_LINK+768_WEDGE+512_STAR+2304_PATH+64_ACTIVE_CYCLE_JETS
STENCILS=1654_CLASSES_HASH_COMMITTED_AND_EXECUTABLY_RECONSTRUCTED;8_LITERAL_FIRST_STENCILS_STORED
DIAGONAL_VALUES=h2:-128,h4:-224/3,h6:-28576/135
ACTIVE_CYCLE_WRITERS=64;H6=-63/8;FIRST_RAW=105/8;SECOND=NONZERO
CTP=DEFINED_Z_PLUS_MINUS_AND_W_MINUS_I_LOG_Z;DIRECT_CONTACT=-SIGMA_ENERGY_HESSIAN_DELTA
RUNTIME_SECONDS=407.049
"""


real_run = subprocess.run
calls = []


def completed_replay(args, *pargs, **kwargs):
    assert len(calls) == 0, "target verifier requested more than one subprocess"
    assert any(str(item).endswith(EXPECTED_DERIVER) for item in args), (
        "target verifier requested an unexpected subprocess")
    calls.append(tuple(str(item) for item in args))
    return subprocess.CompletedProcess(args, 0, stdout=REPLAY_STDOUT, stderr="")


subprocess.run = completed_replay
captured = StringIO()
try:
    with redirect_stdout(captured):
        runpy.run_path(str(TARGET_VERIFIER), run_name="__main__")
finally:
    subprocess.run = real_run

assert len(calls) == 1, "target verifier did not consume completed replay"
output = captured.getvalue()
assert "PASS__GL6CU_V002_PACKET__188/188" in output, output
print(output, end="")
print("PASS__GL6CU_V002_CUSTODY_WITH_COMPLETED_REPLAY__ONE_REPLAY_REUSED")

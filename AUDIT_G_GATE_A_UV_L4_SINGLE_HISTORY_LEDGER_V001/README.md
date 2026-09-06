# Independent hostile audit — L=4 single-history ledger witness

This audit independently reconstructs the one-history calculation without
importing or executing the target verifier.  It attacks four possible errors:

1. a mistaken raw-source coefficient or a source calibration that secretly
   uses the desired ledger result;
2. a circular assignment of the source-owned write term from retained charge;
3. a hidden nonzero transfer across one of the six declared boundary edges;
4. promotion of a controlled F3-MDC witness into a bare-F3 or gravity claim.

The audit can pass only at the stated conditional scope.  It cannot derive
the finite-witness attachment from bare F3, nor does it close the complete
Gate-A owner census.

Run from the repository root:

```text
python3 -B AUDIT_G_GATE_A_UV_L4_SINGLE_HISTORY_LEDGER_V001/independent_reconstruction.py
python3 -B AUDIT_G_GATE_A_UV_L4_SINGLE_HISTORY_LEDGER_V001/verify_audit.py
```

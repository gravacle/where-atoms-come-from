# T-52 GAP ENUMERATION (mechanical, read-only)
Date: 2026-08-21. Method fixed per LANE_T9_AUDIT/METHOD.md; verdicts are for THE REGISTRAR to append.

## AWK USED
awk -F'\t' 'NR==FNR{if(FNR>1)audited[$1]=1;next} FNR>1 && ($4=="FORMAL"||$4=="PROVED"||$4=="MEASURED") && !($1 in audited){print $1"\t"$4}' LANE_T9_AUDIT/T9_carrier_audit.tsv ledger/status_ledger.tsv

## OUTPUT
A-PR	FORMAL
C-90	FORMAL
C-91	FORMAL

## COVERAGE ARITHMETIC
Ledger rows in scope (FORMAL/PROVED/MEASURED): 148. Audit TSV rows: 146.
148 = 145 (audited and still in scope) + 3 (gap above).
One audited row has LEFT the scope since the audit: C-72 is now PARTIAL in the ledger
(its audit row remains in the TSV; status drift, not a gap).

## SIDE FINDING (evidence, not a verdict)
A-PR's GROUNDED cell asserts "no load-bearing citation WITHDRAWN/FAILED/PARTIAL"
(and LANE_T15_ROLES/A_PR_RESTATED.md line 25 records that gate passing), but C-72 —
load-bearing in the EM-world clause — is now PARTIAL. The gate was clean at landing;
it is no longer true of the current ledger.

Per-row carrier citations are in the returned report; gap_rows.tsv beside this file is the raw awk output.

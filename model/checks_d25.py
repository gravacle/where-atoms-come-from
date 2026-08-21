"""FAMILY d25 (T-55): THE D-25 GATE MADE UNBYPASSABLE, AND STAGE (2)'S FIRST REAL ANCHOR GATE.

PROOF_V002 P-DEF-7's finding, made executable: URM.surface refuses undeclared surfaces and is
tested, yet the program's own validators construct RecordSurface directly (model/validate_project.py
lines 4, 20, 44, 53, 54), so no proof number passed through the guard; and D-25 stage (2) — the
source-pinned external anchors carried as GATED checks — was met only by a substring test
("Weller" in provenance).  This block supplies both missing gates:

  1  CONSTRUCTION SCAN (the tool-refusal form of D-25 stage (1)).  An AST scan of the model
     tree FAILS if any file except project_model.py binds, imports, or calls the name
     RecordSurface.  AST, so strings and comments never false-positive and the offense is the
     NAME in code, which also catches aliased imports (RecordSurface as RS).  The scan's zero
     carries two positive controls beside it (D-15): a planted plain offender and a planted
     alias-import offender, both detected in a synthetic tree.
  2  EXTERNAL-ANCHOR GATE (D-25 stage (2), replacing the substring test as its first real
     gate).  The pinned measured azobenzene thermal cis half-life is placed beside the model's
     computed lifetime through URM.surface + ProjectModel.lifetime at the pinned activation
     enthalpy, with the tolerance and semantics stated in each check, and a power control
     showing the gate rejects the historically wrong midpoint convention (the C-69 demotion
     note's cautionary tale).

CLAIM ROWS: D-25 (DOCTRINE, enforcement stages 1-2); the C-69 demotion note (what the lifetime
law's agreement does and does not establish); P-DEF-7 (the bypass this closes).
SEALED SOURCES: LANE_T41_EXTERNAL/CITATIONS.md (source-pinned anchors, sealed);
ledger/status_ledger.tsv rows D-25 and C-69; PROOF_V002.md P-DEF-7.
OWNERS: the measured datum and the pinned Eyring parameters are the literature's (Thermal
Half-Lives of Azobenzene Derivatives, ACS Cent. Sci., PMC9951306); the Eyring/Arrhenius
formalism is standard and borrowed; the scan instrument, the gate semantics, and the power
control are the program's.

PINNED EXTERNAL ANCHORS (D-8: these literals ARE the anchors, stated as such; every model-side
number on a decision path below is computed through the URM, never written down):
  t1/2 = 1.2e5 s        measured thermal cis->trans half-life, 1.4 d in benzene at 35 C
  dH   = 0.915 eV       pinned activation enthalpy, 21.1 kcal/mol
  dS   = -50.2 J/mol/K  pinned activation entropy
  (all three: ACS Cent. Sci. PMC9951306 via LANE_T41_EXTERNAL/CITATIONS.md item 1)
  dE   = 0.60 eV        cis-trans energy difference, the corpus's azobenzene value
                        (model/validate_project.py T-33 row; enters only the negligible
                        reverse rate)

Importable standalone; runtime well under one second."""
import ast
import math
import os
import shutil
import sys
import tempfile

_HERE = os.path.dirname(os.path.abspath(__file__))
if _HERE not in sys.path:
    sys.path.insert(0, _HERE)

import grounded as G

# ---- pinned external anchors (see header) ------------------------------------------------
AZO_T_HALF_S = 1.2e5          # measured t1/2, benzene, 35 C
AZO_DH_EV = 0.915             # pinned activation enthalpy
AZO_DS_J_MOL_K = -50.2        # pinned activation entropy
AZO_T_K = 308.15              # 35 C, the measurement's own condition
AZO_DE_EV = 0.60              # corpus cis-trans energy difference
# ---- defined SI constants (exact, CODATA/SI; not anchors) --------------------------------
EV = 1.602176634e-19          # J per eV, exact
H_PLANCK = 6.62607015e-34     # J s, exact
R_GAS = 8.314462618           # J/mol/K (= KB * N_A, exact)
LN2 = math.log(2.0)

_TARGET_NAME = "RecordSurface"


# ==========================================================================================
# 1  THE CONSTRUCTION SCAN — D-25 stage (1) in tool-refusal form
# ==========================================================================================
def scan_direct_constructions(root=None, allowed=("project_model.py",)):
    """AST-scan every .py under root (default: this file's directory, the model tree) for the
    name RecordSurface appearing in CODE — an import (any alias), a bare-name reference, or an
    attribute access — outside the allowed basenames.  Returns [(relpath, lineno, kind), ...],
    empty when the tree is clean.  A file that does not parse is itself an offense: what
    cannot be audited is refused."""
    root = root or _HERE
    offenders = []
    for dirpath, dirnames, filenames in os.walk(root):
        dirnames[:] = sorted(d for d in dirnames if d != "__pycache__")
        for fn in sorted(filenames):
            if not fn.endswith(".py") or fn in allowed:
                continue
            path = os.path.join(dirpath, fn)
            with open(path, encoding="utf-8") as fh:
                src = fh.read()
            try:
                tree = ast.parse(src)
            except SyntaxError as e:
                offenders.append((os.path.relpath(path, root), int(e.lineno or 0),
                                  "unparseable"))
                continue
            for node in ast.walk(tree):
                kind = None
                if isinstance(node, ast.ImportFrom):
                    if any(a.name == _TARGET_NAME for a in node.names):
                        kind = "import"
                elif isinstance(node, ast.Name) and node.id == _TARGET_NAME:
                    kind = "reference"
                elif isinstance(node, ast.Attribute) and node.attr == _TARGET_NAME:
                    kind = "attribute"
                if kind:
                    offenders.append((os.path.relpath(path, root),
                                      int(getattr(node, "lineno", 0)), kind))
    return offenders


def _planted_tree(files):
    """A throwaway directory holding the given {basename: source} files, for the scan's
    positive controls.  Caller removes it."""
    d = tempfile.mkdtemp(prefix="d25_scan_control_")
    for fn, src in files.items():
        with open(os.path.join(d, fn), "w", encoding="utf-8") as fh:
            fh.write(src)
    return d


# ==========================================================================================
# 2  THE CHECK BLOCK — validate_geometry.py idiom
# ==========================================================================================
def run_d25_checks(check, root=None):
    """All family-d25 gates through the caller's check(name, cond, detail).  root overrides
    the scanned tree (default: this file's directory) so the block can be pointed at a
    candidate tree; the chained validator uses the default."""
    from project_model import ProjectModel, URM

    M = ProjectModel()

    # ---------------------------------- D-25 stage (1): the scan and its controls (D-15)
    off = scan_direct_constructions(root=root)
    check("D-25 scan: no RecordSurface binding outside project_model.py in the model tree",
          off == [],
          "clean" if off == [] else "OFFENDERS " + "; ".join(
              f"{p}:{ln} ({k})" for p, ln, k in off))

    plain = _planted_tree({
        "offender.py": "from project_model import " + _TARGET_NAME + "\n"
                       "s = " + _TARGET_NAME + "('x', 'y', 1.0, 2.0, 300.0, 1e9)\n",
        "clean.py": "from project_model import URM\n",
    })
    try:
        got = scan_direct_constructions(root=plain)
    finally:
        shutil.rmtree(plain, ignore_errors=True)
    check("D-25 scan control: planted direct construction is detected (positive control "
          "beside the zero)",
          any(p == "offender.py" for p, _, _ in got)
          and all(p != "clean.py" for p, _, _ in got),
          f"flagged {sorted(set(p for p, _, _ in got))}, clean.py unflagged")

    alias = _planted_tree({
        "sly.py": "from project_model import " + _TARGET_NAME + " as RS\n"
                  "s = RS('x', 'y', 1.0, 2.0, 300.0, 1e9)\n",
    })
    try:
        got2 = scan_direct_constructions(root=alias)
    finally:
        shutil.rmtree(alias, ignore_errors=True)
    check("D-25 scan probe: aliased import (as RS) is detected — the scan keys on the "
          "imported name, not the call text (definition, no shortcut)",
          any(p == "sly.py" and k == "import" for p, _, k in got2),
          f"flagged {got2}")

    # ---------------------------------- the refusal beyond the gated cases (API fidelity)
    refused_ws = False
    try:
        URM.surface("mystery device", "unknown", 1e-20, 1e-19, 300.0, 1e9, provenance="   ")
    except ValueError:
        refused_ws = True
    check("D-25 refusal probe: whitespace-only provenance is refused (beyond the "
          "validator's gated blank case)", refused_ws, "ValueError raised")

    refused_lc = False
    try:
        URM.surface("toy torus", "stabiliser", 0.0, 0.0, 0.0, 0.0, tier="corner",
                    provenance="def-a")
    except ValueError:
        refused_lc = True
    check("D-25 refusal probe: corner tier refuses lowercase 'def-a' — the DEF-A "
          "self-declaration is exact", refused_lc, "ValueError raised")

    okw = URM.surface("Azobenzene", "photoisomerisation", AZO_DE_EV * EV, AZO_DH_EV * EV,
                      AZO_T_K, 1e12)
    okc = URM.surface("toy torus", "stabiliser", 0.0, 0.0, 0.0, 0.0, tier="corner",
                      provenance="DEF-A")
    check("D-25 positive control beside the refusals: declared surfaces construct through "
          "the gate (world via registry, corner via DEF-A)",
          okw.provenance is not None and "PMC9951306" in okw.provenance
          and okc.tier == "corner" and okc.provenance == "DEF-A",
          f"world provenance pinned: {okw.provenance[:52]}...; corner tier {okc.tier}")

    # ---------------------------------- D-25 stage (2): the first real external-anchor gate
    # Every model lifetime below is computed through URM.surface + ProjectModel.lifetime at
    # the pinned dH and the measurement's own temperature; only f0 varies.
    def azo_tau(E_b_eV, f0, note):
        s = URM.surface("Azobenzene", "photoisomerisation", AZO_DE_EV * EV, E_b_eV * EV,
                        AZO_T_K, f0, provenance=note)
        return M.lifetime(s)

    pin = ("pinned: ACS Cent. Sci. PMC9951306 via LANE_T41_EXTERNAL/CITATIONS.md item 1; "
           "dH = 0.915 eV, T = 308.15 K")
    tau_meas = AZO_T_HALF_S / LN2   # first-order kinetics: mean lifetime = t1/2 / ln 2

    tau_hi = azo_tau(AZO_DH_EV, 1e9, pin + "; envelope edge f0 = 1e9 /s")
    tau_lo = azo_tau(AZO_DH_EV, 1e13, pin + "; envelope edge f0 = 1e13 /s")
    check("D-25 anchor (envelope): measured azobenzene mean lifetime lies inside the model "
          "envelope at pinned dH = 0.915 eV over f0 in [1e9, 1e13] /s",
          tau_lo <= tau_meas <= tau_hi,
          f"tau_meas = t1/2/ln2 = {tau_meas:.3e} s in [{tau_lo:.3e}, {tau_hi:.3e}] s. "
          "SEMANTICS: containment in a four-decade attempt-frequency envelope — weak "
          "contact; consistency of the pinned barrier with the measured lifetime for some "
          "physical prefactor, no prediction. The f0 class spans condensed-phase "
          "prefactors down-shifted by the pinned strongly negative dS (Eyring gives "
          "~1.5e10 /s, inside the class).")

    # the Eyring point: f0 computed from the SAME pinned source's activation entropy —
    # computed here, never written down
    f0_eff = (G.KB * AZO_T_K / H_PLANCK) * math.exp(AZO_DS_J_MOL_K / R_GAS)
    tau_model = azo_tau(AZO_DH_EV, f0_eff,
                        pin + f"; Eyring prefactor from pinned dS = {AZO_DS_J_MOL_K} J/mol/K")
    t_half_model = tau_model * LN2
    dec = math.log10(t_half_model / AZO_T_HALF_S)
    check("D-25 anchor (Eyring point): model t1/2 at the pinned (dH, dS) within ONE DECADE "
          "of the measured t1/2 — the stated tolerance",
          abs(dec) <= 1.0,
          f"model t1/2 = {t_half_model:.3e} s vs measured {AZO_T_HALF_S:.1e} s: "
          f"{dec:+.2f} decades (factor {AZO_T_HALF_S/t_half_model:.1f} — the pinned "
          "source's own internal spread between its Eyring fit and its quoted half-life). "
          "SEMANTICS: order-of-magnitude consistency when the pinned source supplies BOTH "
          "parameters. It does NOT establish prediction — the prefactor comes from the same "
          "source — and per the C-69 demotion note the closed-form agreement is an identity "
          "of the construction (the rates ARE the inserted jump strengths); the model's own "
          "content here is that lifetime() reads the record's Liouvillian mode.")

    # power control (D-15): the historically wrong midpoint convention must FAIL this gate.
    # The old convention silently reduced the escape barrier by dE/2 — the 5-6 order
    # azobenzene contradiction that forced the correction (C-69 demotion note).
    tau_mid = azo_tau(AZO_DH_EV - AZO_DE_EV / 2.0, f0_eff,
                      "CONTROL: pre-correction midpoint convention, barrier reduced by "
                      "dE/2 — historically wrong, must fail the decade gate (C-69)")
    dec_mid = math.log10(tau_mid * LN2 / AZO_T_HALF_S)
    check("D-25 anchor power control: the pre-correction midpoint convention FAILS the same "
          "decade gate — the gate can reject a wrong convention",
          abs(dec_mid) > 1.0,
          f"midpoint-convention t1/2 = {tau_mid * LN2:.3e} s: {dec_mid:+.2f} decades off "
          "(the 5-6 order contradiction the C-69 demotion note records)")

    # D-8 cross-check: the Liouvillian-mode lifetime on the decision path equals the closed
    # form computed independently here — the anchor gate rests on the mode extraction.
    kT = G.KB * AZO_T_K
    gu = f0_eff * math.exp(-(AZO_DH_EV * EV) / kT)
    gd = f0_eff * math.exp(-((AZO_DH_EV + AZO_DE_EV) * EV) / kT)
    rel = abs(tau_model - 1.0 / (gu + gd)) / (1.0 / (gu + gd))
    check("D-25 anchor cross-check: the gated lifetime is the record's own Liouvillian mode "
          "(equals the independent closed form to <1e-9 relative)",
          rel < 1e-9, f"rel err {rel:.2e}")


# ==========================================================================================
# standalone runner
# ==========================================================================================
if __name__ == "__main__":
    n_pass = 0
    n_fail = 0

    def check(name, cond, detail=""):
        global n_pass, n_fail
        if cond:
            n_pass += 1
            print(f"  PASS  {name}  {detail}")
        else:
            n_fail += 1
            print(f"  FAIL  {name}  {detail}")

    print("VALIDATE FAMILY d25 (T-55): THE D-25 GATE MADE UNBYPASSABLE")
    print("=" * 78)
    root = sys.argv[1] if len(sys.argv) > 1 else None
    run_d25_checks(check, root=root)
    print("=" * 78)
    print(f"  D25: {n_pass} PASS, {n_fail} FAIL")
    sys.exit(0 if n_fail == 0 else 1)

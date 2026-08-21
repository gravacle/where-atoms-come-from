#!/usr/bin/env python3
"""THE PROOF GATE — read PROOF_V002.md against the RECORD and refuse it when the two disagree.

  ./replicate/check_proof.py                     check PROOF_V002.md at the repo root
  ./replicate/check_proof.py path/to/draft.md    check some other draft

T-16's DONE_WHEN is "PROOF_V002.md, every step citing a ledger row, no step resting on a single
carrier unless marked". This file is what makes that impossible to violate SILENTLY. The program's
standing rule, from the principal: A GUARD MUST BE A TOOL REFUSAL, NEVER A CHECKBOX. Every guard in
this program that was written as a checklist item made the problem LOOK handled and changed nothing
— H-3 sat open for the life of the program while 162 rows stood marked PROVED about a stipulation.

Nothing here is compared against a number or a status typed into this file (D-8). The document is
compared against `ledger/status_ledger.tsv` and `LANE_T9_AUDIT/T9_carrier_audit.tsv`, which are the
record. When those cannot be read the tool exits NONZERO instead of passing quietly: a gate that
reports PASS because it could not find the ledger is worse than no gate at all.

THE RULES, and why each one exists:

  R1   every "### P-" block parses and carries all five fields. A block missing its gate or its
       grounding is a claim with no way to be wrong.
  R2   every cited ID is a row in the ledger. A citation to a row that does not exist is prose
       wearing a citation's clothes.
  R3   no cited row is WITHDRAWN, FAILED or RECLASSIFIED. The register is append-only precisely so
       that a retracted claim stays visible — and stays uncitable.
  R4   the status PRINTED in the proof equals the LEDGER's status. This is the failure mode the
       tool exists for: the ledger moves, the proof does not, and the document keeps asserting a
       status that was withdrawn under it.
  R5   the carrier mark printed equals the T-9 audit's verdict, or is UNAUDITED when the audit
       never reached that row. Rows registered after the audit are UNAUDITED, which is NOT a mark.
  R6   a block that rests on no TWO-CARRIER row must SAY SO in scope. T-9's whole finding is that
       most of this program stands on one carrier; a proof that does not carry that forward is
       claiming more than the record holds.
  R7   the model cell names a function that exists under model/. The architecture the principal
       fixed for T-16 is that THE MODEL IS THE PROOF; a step with no model function behind it is
       narration and must say `none`.
  R8   the gate cell names a check that actually fires. A gate nobody runs is a checkbox.
  R9   a block whose prose says PROVED cites a row the ledger calls PROVED. FORMAL is real
       mathematics about our own stipulated definition and says nothing about the world (H-3).
  R10  D-1, the banned comparisons. Classical gravity is not expected at the record level and is
       never the test. Newton may be named in ONE section, and only there.
  R11  the summary line, printed on every run, PASS or FAIL, because the true state of the document
       must be impossible not to see.

Exit: 0 clean · 1 a rule fired · 2 this tool's own assumptions are broken.
"""
import ast
import csv
import os
import re
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
LEDGER  = os.path.join(ROOT, 'ledger', 'status_ledger.tsv')
PLAN    = os.path.join(ROOT, 'ledger', 'plan.tsv')
VOCAB   = os.path.join(ROOT, 'ledger', 'STATUS_VOCAB.tsv')
AUDIT   = os.path.join(ROOT, 'LANE_T9_AUDIT', 'T9_carrier_audit.tsv')
ONFIND  = os.path.join(ROOT, 'replicate', 'on_finding.sh')
MODELDIR = os.path.join(ROOT, 'model')
DEFAULT_PROOF = os.path.join(ROOT, 'PROOF_V002.md')

# THE FORMAT CONTRACT. A claim block is a "### P-" heading, prose, and a two-column pipe table
# carrying exactly these five fields. The order of the tuple is the order they are reported in.
FIELDS = ('model', 'gate', 'grounding', 'rows', 'scope')
NONE = 'none'                       # the literal word a block writes when there is no model / no gate

# A row in one of these three states was tested and did not hold, or was retracted. It stays in the
# ledger forever (rows are never deleted) and it may never be leaned on again.
DEAD = ('WITHDRAWN', 'FAILED', 'RECLASSIFIED')

# The one status that means a test any physicist anywhere can run against their own data. Both this
# and DEAD are checked against ledger/STATUS_VOCAB.tsv at startup: if the closed vocabulary moves,
# this tool must be fixed, not left quietly agreeing with a word the ledger no longer uses.
PROVED = 'PROVED'

# The carrier marks. UNAUDITED is not a verdict the T-9 audit ever writes — it is what a row that
# the audit never reached carries, and it is NOT a two-carrier mark.
TWO_CARRIER = 'TWO-CARRIER'
UNAUDITED = 'UNAUDITED'
MARKS = (TWO_CARRIER, 'SINGLE-CARRIER', 'NOT-CARRIER-SHAPED', UNAUDITED)

# The prefix a block must carry in `scope` when nothing it cites has a second carrier. Written with
# an explicit escape so the em dash cannot be silently replaced by a hyphen in an editor.
SINGLE_PREFIX = 'SINGLE-CARRIER — '

MIDDOT = '·'                   # the separator inside the rows cell

RULES = [
    ('R1',  'every P- block parses and carries all five fields'),
    ('R2',  'every cited ID is a row in the ledger'),
    ('R3',  'no cited row is WITHDRAWN, FAILED or RECLASSIFIED'),
    ('R4',  'the status printed in the proof is the ledger\'s status'),
    ('R5',  'the carrier mark printed is the T-9 audit\'s verdict'),
    ('R6',  'a block with no TWO-CARRIER row says SINGLE-CARRIER in scope'),
    ('R7',  'every model cell names a function that exists under model/'),
    ('R8',  'every gate cell names a check that actually fires'),
    ('R9',  'PROVED in the prose is backed by a PROVED row'),
    ('R10', 'no banned comparison outside the comparison section (D-1)'),
]


def fatal(msg):
    """Die loudly. Reserved for this tool's OWN broken assumptions — a missing ledger, an audit
       file with a verdict we do not know how to read — never for a finding about the document."""
    sys.stderr.write('CHECK_PROOF REFUSES TO RUN\n\n  ' + msg.replace('\n', '\n  ') + '\n')
    sys.exit(2)


# ----------------------------------------------------------------------------- the record

def read_tsv(path, what):
    """QUOTE_NONE, deliberately: these files are TSVs and a stray double quote inside an evidence
       cell must stay a character, not turn into a quoting rule that swallows the next row."""
    if not os.path.exists(path):
        fatal('CANNOT FIND %s\n    %s\nThe proof is checked AGAINST the record. Without the record\n'
              'there is nothing to check against, and this tool will not report PASS.' % (what, path))
    try:
        with open(path, newline='', encoding='utf-8') as f:
            rows = list(csv.reader(f, delimiter='\t', quoting=csv.QUOTE_NONE))
    except Exception as exc:                                    # unreadable is not the same as clean
        fatal('CANNOT READ %s\n    %s\n%s' % (what, path, exc))
    rows = [r for r in rows if r and any(c.strip() for c in r)]
    if len(rows) < 2:
        fatal('%s has a header and no rows: %s' % (what, path))
    return rows[0], rows[1:]


def load_ledger():
    """ID -> STATUS, by COLUMN NAME. Reading column 4 by position is how a tool starts checking the
       wrong field the day a column is inserted."""
    hdr, rows = read_tsv(LEDGER, 'the status ledger')
    try:
        i_id, i_st = hdr.index('ID'), hdr.index('STATUS')
    except ValueError:
        fatal('the status ledger has no ID/STATUS columns; its header reads: %s' % hdr)
    out = {}
    for r in rows:
        if len(r) > max(i_id, i_st) and r[i_id].strip():
            out[r[i_id].strip()] = r[i_st].strip()
    return out


def load_plan():
    """The plan is the second register: a standing-debt block rests on the TASK that owes the work,
       and `T-38 (TODO, PLAN)` is a citation with exactly as much force as a ledger row. A plan task
       is never carrier evidence — PLAN is a tier, not a mark — so R6 still bites on such a block."""
    hdr, rows = read_tsv(PLAN, 'the plan')
    try:
        i_id, i_st = hdr.index('ID'), hdr.index('STATUS')
    except ValueError:
        fatal('the plan has no ID/STATUS columns; its header reads: %s' % hdr)
    out = {}
    for r in rows:
        if len(r) > max(i_id, i_st) and r[i_id].strip():
            out[r[i_id].strip()] = r[i_st].strip()
    if not out:
        fatal('ledger/plan.tsv holds no task. R2 would then refuse every debt citation.')
    return out


def load_vocab():
    """The closed status vocabulary, used to confirm THIS TOOL still speaks the ledger's language."""
    _, rows = read_tsv(VOCAB, 'the closed status vocabulary')
    return {r[0].strip() for r in rows if r[0].strip()}


def load_audit():
    """ID -> the T-9 carrier verdict. Column 1 is the ID, column 2 the verdict, as sealed."""
    _, rows = read_tsv(AUDIT, 'the T-9 carrier audit')
    out = {}
    for r in rows:
        if r[0].strip():
            out[r[0].strip()] = (r[1].strip() if len(r) > 1 else '')
    unknown = sorted(set(out.values()) - set(MARKS))
    if unknown:
        fatal('the T-9 audit carries a verdict this tool does not know: %s\n'
              'R5 was written against a vocabulary that has since moved. FIX THIS TOOL. Do not\n'
              'let it keep passing documents it can no longer read.' % unknown)
    return out


def load_banned_pattern():
    """D-1's pattern is READ OUT OF replicate/on_finding.sh, never re-typed here. Two copies of a
       ban drift apart, and the copy that drifts is always the one that stops firing.

       on_finding.sh then filters its hits with a `grep -vi` list, because the LEDGER legitimately
       contains rows that DESCRIBE the ban. That lexical escape is deliberately NOT reused here: in
       a proof the exemption is STRUCTURAL — the one comparison section — so a line that names a
       banned comparison anywhere else is a finding, whatever other words it happens to contain."""
    if not os.path.exists(ONFIND):
        fatal('CANNOT FIND replicate/on_finding.sh (%s); D-1\'s pattern lives there and this tool\n'
              'reads it rather than keeping a second copy.' % ONFIND)
    src = open(ONFIND, encoding='utf-8').read()
    m = re.search(r'^PAT="(.+)"\s*$', src, re.M)
    if not m:
        fatal('replicate/on_finding.sh no longer defines PAT="..." on one line.\n'
              'D-1\'s scan is single-sourced from that file. FIX ONE OF THE TWO, do not fork it.')
    try:
        return re.compile(m.group(1), re.I)
    except re.error as exc:
        fatal('replicate/on_finding.sh\'s PAT is not a pattern python can compile: %s' % exc)


def load_model_defs():
    """Every `def <name>` under model/. The proof's model cell must land on one of these: T-16's
       architecture is that THE MODEL IS THE PROOF, so a claim whose carrier function does not
       exist is narration, and narration is required to say `none`."""
    if not os.path.isdir(MODELDIR):
        fatal('CANNOT FIND model/ (%s). The proof narrates the model; without it R7 is blind.' % MODELDIR)
    names = set()
    pat = re.compile(r'^[ \t]*def[ \t]+([A-Za-z_]\w*)[ \t]*\(', re.M)
    for base, _, files in os.walk(MODELDIR):
        for fn in files:
            if fn.endswith('.py'):
                try:
                    names.update(pat.findall(open(os.path.join(base, fn), encoding='utf-8').read()))
                except Exception as exc:
                    fatal('CANNOT READ %s: %s' % (os.path.join(base, fn), exc))
    if not names:
        fatal('no `def` found anywhere under model/. That cannot be right; R7 would pass nothing.')
    return names


def check_names(path):
    """The gate names a block may cite are exactly the strings a validator hands to check(). They
       are read from the SOURCE by parsing it. Three names carry them across the validators
       — check(), and chk()/chkb() in model/validate_formation.py — and all three are read — a list of names typed into this file would go stale
       the first time a validator was edited, which is precisely the failure R4 exists to catch.

       Returns (exact names, regexes for f-string names). An f-string gate name is not literally
       present in the source, so it is matched by its template rather than pretended not to exist."""
    try:
        tree = ast.parse(open(path, encoding='utf-8').read(), filename=path)
    except SyntaxError as exc:
        fatal('%s does not parse as python, so its check names cannot be read: %s' % (path, exc))
    exact, templates = set(), []
    for node in ast.walk(tree):
        if not isinstance(node, ast.Call) or not node.args:
            continue
        fn = node.func
        name = getattr(fn, 'id', None) or getattr(fn, 'attr', None)
        if name not in ('check', 'chk', 'chkb'):
            continue
        first = node.args[0]
        if isinstance(first, ast.Constant) and isinstance(first.value, str):
            exact.add(first.value)
        elif isinstance(first, ast.JoinedStr):
            parts = [re.escape(v.value) if isinstance(v, ast.Constant) else '.+'
                     for v in first.values]
            templates.append(re.compile('^' + ''.join(parts) + '$'))
    return exact, templates


# ----------------------------------------------------------------------------- the document

HEADING = re.compile(r'^(#{1,6})\s+(.*)$')
BLOCK_HEAD = re.compile(r'^###\s+(P-\S+)')
# A SECTION heading, never a claim block's. Blocks are `### P-`, and a block whose headline happens
# to contain the word must not silently open a second exemption — the exemption is structural.
COMPARISON_HEAD = re.compile(r'^##(?!#)\s+.*COMPARISON', re.I)
CITE = re.compile(r'^([A-Za-z][A-Za-z0-9.]*-[A-Za-z0-9.\-]+)\s*(?:\((.*)\))?$')


class Block(object):
    def __init__(self, bid, line):
        self.id = bid
        self.line = line            # 1-indexed, so a violation can be walked to
        self.prose = []
        self.fields = {}            # name -> cell text, only from two-column rows
        self.cites = []             # (raw, id, printed status or None, printed mark or None)


def split_cells(line):
    """A pipe row is a two-column row only if it yields exactly two cells. A row that does not is
       simply not one of the contract's field rows, and the field it meant to carry is reported
       MISSING by R1 — which is the honest reading, and refuses a three-column table by itself."""
    s = line.strip()
    if not s.startswith('|'):
        return None
    # A cell may legitimately CONTAIN a pipe -- |lambda| is an operator norm, not a column break --
    # and markdown escapes it as \|. Honour the escape, or the field vanishes and R1 complains
    # about the wrong thing.
    parts, buf, i = [], [], 0
    while i < len(s):
        ch = s[i]
        if ch == '\\' and i + 1 < len(s) and s[i + 1] == '|':
            buf.append('|'); i += 2; continue
        if ch == '|':
            parts.append(''.join(buf)); buf = []; i += 1; continue
        buf.append(ch); i += 1
    parts.append(''.join(buf))
    if parts and not parts[0].strip():
        parts = parts[1:]
    if parts and not parts[-1].strip():
        parts = parts[:-1]
    cells = [c.strip() for c in parts]
    return cells if len(cells) == 2 else None


def parse_blocks(lines):
    blocks, cur = [], None
    for n, raw in enumerate(lines, 1):
        head = BLOCK_HEAD.match(raw)
        if head:
            cur = Block(head.group(1), n)
            blocks.append(cur)
            continue
        if cur is None:
            continue
        if re.match(r'^#{1,3}\s', raw):                  # any heading at ### or above closes a block
            cur = None
            continue
        cells = split_cells(raw)
        if cells is None:
            if raw.strip():
                cur.prose.append(raw.strip())
            continue
        key = cells[0].strip().strip('*').strip('`').strip().lower()
        if not key or set(key) <= set('-: '):            # the empty header row and the separator
            continue
        if key in FIELDS and key not in cur.fields:
            cur.fields[key] = cells[1].strip()
    return blocks


def parse_cites(cell):
    """`C-71 (PROVED, TWO-CARRIER) · C-72 (PROVED, TWO-CARRIER)` — the form the registrar leaves
       behind once the statuses and marks have been filled in mechanically. Backticks around an ID
       are cosmetic and are stripped wherever they fall: a formatting choice may not defeat a gate."""
    out = []
    for part in cell.split(MIDDOT):
        # backticks are a formatting choice; they must not be a way to slip past the parser
        s = part.replace('`', ' ').strip()
        if not s:
            continue
        m = CITE.match(s)
        if not m:
            out.append((s, None, None, None))
            continue
        rid, inner = m.group(1), m.group(2)
        if inner is None:
            out.append((s, rid, None, None))
        elif ',' in inner:
            st, mk = inner.split(',', 1)
            out.append((s, rid, st.strip(), mk.strip()))
        else:
            out.append((s, rid, inner.strip(), None))
    return out


def comparison_exempt(lines):
    """The ONE place Newton may be named. The exemption runs from the comparison heading to the
       next heading at the same level or above — a section, not a licence for the rest of the file.
       Returns (set of exempt line numbers, list of comparison heading line numbers)."""
    exempt, heads = set(), []
    for n, raw in enumerate(lines, 1):
        if COMPARISON_HEAD.match(raw):
            heads.append(n)
    for start in heads:
        level = len(re.match(r'^#+', lines[start - 1]).group(0))
        for n in range(start, len(lines) + 1):
            if n > start:
                h = HEADING.match(lines[n - 1])
                if h and len(h.group(1)) <= level:
                    break
            exempt.add(n)
    return exempt, heads


# ----------------------------------------------------------------------------- the gate

PLAN_TIER = 'PLAN'


def main():
    args = [a for a in sys.argv[1:] if not a.startswith('-')]
    if any(a in ('-h', '--help') for a in sys.argv[1:]):
        sys.stdout.write(__doc__)
        return 0
    proof = os.path.abspath(args[0]) if args else DEFAULT_PROOF
    if not os.path.exists(proof):
        fatal('CANNOT FIND the proof: %s\nT-16 is not done until this file exists.' % proof)

    ledger = load_ledger()
    planned = load_plan()
    vocab = load_vocab()
    missing = [s for s in DEAD + (PROVED,) if s not in vocab]
    if missing:
        fatal('the closed vocabulary no longer contains %s.\n'
              'R3 and R9 were written against a vocabulary that has since moved. FIX THIS TOOL.'
              % missing)
    audit = load_audit()
    banned = load_banned_pattern()
    defs = load_model_defs()

    text = open(proof, encoding='utf-8').read()
    lines = text.split('\n')
    blocks = parse_blocks(lines)

    fired = []                       # (rule, block id or '-', message)
    def hit(rule, where, msg):
        fired.append((rule, where, msg))

    if not blocks:
        hit('R1', '-', 'no "### P-" claim block anywhere in the document — under the V002 contract '
                       'this is not a proof, it is prose')

    for b in blocks:
        # R1 — the five fields. A block missing its gate or its grounding is a claim with no way
        # to be wrong, which is the one thing a proof step may not be.
        for f in FIELDS:
            if f not in b.fields:
                hit('R1', b.id, 'line %d: no **%s** row in its table' % (b.line, f))
            elif not b.fields[f]:
                hit('R1', b.id, 'line %d: **%s** is empty' % (b.line, f))

        b.cites = parse_cites(b.fields.get('rows', ''))
        if 'rows' in b.fields and b.fields['rows'] and not b.cites:
            hit('R1', b.id, 'line %d: **rows** carries no citation' % b.line)

        cited_ids = []
        for raw, rid, printed_st, printed_mk in b.cites:
            if rid is None:
                hit('R2', b.id, 'cannot read the citation %r — expected  ID (STATUS, MARK)' % raw)
                continue
            cited_ids.append(rid)

            # R2 — the ledger is the register of what this program claims. A citation to something
            # that is not in it is not a citation.
            if rid not in ledger:
                if rid in planned:
                    # a standing debt cites the task that owes it. Same discipline, other register.
                    if printed_st is None:
                        hit('R4', b.id, '%s is cited with no (STATUS, TIER)' % rid)
                    elif printed_st != planned[rid]:
                        hit('R4', b.id, '%s is printed %s; the plan says %s'
                            % (rid, printed_st, planned[rid]))
                    if printed_mk != PLAN_TIER:
                        hit('R5', b.id, '%s is a plan task and must be printed (%s, %s)'
                            % (rid, planned[rid], PLAN_TIER))
                    continue
                hit('R2', b.id, '%s is in neither ledger/status_ledger.tsv nor ledger/plan.tsv' % rid)
                continue

            true_st = ledger[rid]
            true_mk = audit.get(rid, UNAUDITED)

            # R3 — rows are appended and never deleted so a retraction stays VISIBLE. Visible, and
            # uncitable.
            if true_st in DEAD:
                hit('R3', b.id, '%s is %s in the ledger and may not be leaned on' % (rid, true_st))

            # R4 — THE REASON THIS TOOL EXISTS. The ledger moves; a document does not move with it.
            if printed_st is None:
                hit('R4', b.id, '%s is cited with no (STATUS, MARK) — the registrar\'s fill is '
                                'missing, so the printed status cannot be checked at all' % rid)
            elif printed_st != true_st:
                hit('R4', b.id, '%s is printed %s; the ledger says %s' % (rid, printed_st, true_st))

            # R5 — a mark is a verdict of the T-9 audit or it is UNAUDITED. Nothing else.
            if printed_mk is None:
                if printed_st is not None:
                    hit('R5', b.id, '%s is printed with no carrier mark' % rid)
            elif printed_mk not in MARKS:
                hit('R5', b.id, '%s is printed %s, which is not a carrier mark' % (rid, printed_mk))
            elif printed_mk != true_mk:
                where = ('the T-9 audit says %s' % true_mk if rid in audit else
                         'the T-9 audit never reached %s, so it is %s' % (rid, UNAUDITED))
                hit('R5', b.id, '%s is printed %s; %s' % (rid, printed_mk, where))

        # R6 — T-9's finding is that most of this program stands on ONE carrier. A block resting on
        # nothing better must carry that forward wherever it is quoted, or the quote overclaims.
        has_two = any(audit.get(r) == TWO_CARRIER for r in cited_ids)
        scope = b.fields.get('scope', '')
        if not has_two and not scope.startswith(SINGLE_PREFIX):
            hit('R6', b.id, 'no cited row carries %s, so **scope** must begin %r' %
                (TWO_CARRIER, SINGLE_PREFIX))

        # R7 — the model cell. THE MODEL IS THE PROOF (the principal, T-16 architecture).
        model = b.fields.get('model', '').replace('`', ' ').strip()
        if model and model.lower() != NONE:
            leaf = model.rstrip('()').split('.')[-1].strip()
            if leaf not in defs:
                hit('R7', b.id, '**model** %s: no `def %s` anywhere under model/ — if this claim is '
                                'narration the cell must read `%s`' % (model, leaf, NONE))

        # R8 — the gate cell. A gate that does not fire is a checkbox, and this program has paid
        # for that lesson more than once.
        # backticks are cosmetic here too; the gate is checked, not the markdown around it
        gate = b.fields.get('gate', '').replace('`', ' ').strip()
        gate = re.sub(r'\s*::\s*', ' :: ', gate)
        if gate and gate.lower() != NONE:
            if '::' not in gate:
                hit('R8', b.id, '**gate** %r is neither `%s` nor `<file> :: <check name>`' % (gate, NONE))
            else:
                gpath, gname = [x.strip() for x in gate.split('::', 1)]
                full = gpath if os.path.isabs(gpath) else os.path.join(ROOT, gpath)
                if not os.path.exists(full):
                    hit('R8', b.id, '**gate** names %s, which does not exist' % gpath)
                else:
                    exact, templates = check_names(full)
                    if gname not in exact and not any(t.match(gname) for t in templates):
                        hit('R8', b.id, '**gate** %s :: %r is not a check() in that file' % (gpath, gname))

        # R9 — FORMAL is real mathematics about OUR OWN stipulated definition and says nothing about
        # the world by itself (H-3). The word PROVED in prose has to be earned in the ledger.
        prose = ' '.join(b.prose)
        if re.search(r'\bPROVED\b', prose):
            if not any(ledger.get(r) == PROVED for r in cited_ids):
                hit('R9', b.id, 'the prose says PROVED and no cited row is PROVED in the ledger')

    # R10 — D-1. Classical gravity is not expected at the record level and is never the test. The
    # pattern comes from replicate/on_finding.sh; the only exemption is the comparison section.
    exempt, heads = comparison_exempt(lines)
    if len(heads) > 1:
        hit('R10', '-', 'more than one COMPARISON section (lines %s) — the exemption is ONE section, '
                        'not a licence spread through the document' %
                        ', '.join(str(h) for h in heads))
    for n, raw in enumerate(lines, 1):
        if n in exempt:
            continue
        m = banned.search(raw)
        if m:
            hit('R10', '-', 'line %d: banned comparison %r in: %s' % (n, m.group(0), raw.strip()[:100]))

    # ------------------------------------------------------------------- the report
    rel = os.path.relpath(proof, ROOT) if proof.startswith(ROOT) else proof
    print('PROOF GATE — %s' % rel)
    print('=' * 66)
    struck = {r for r, _, _ in fired}
    for rule, desc in RULES:
        if rule in struck:
            for r, where, msg in fired:
                if r == rule:
                    print('  FAIL  %-4s %-10s %s' % (r, where, msg))
        else:
            print('  ok    %-4s %s' % (rule, desc))

    # R11 — printed on EVERY run, pass or fail. The grounding debt taught this program that a state
    # you have to go and look up is a state nobody looks up.
    distinct = sorted({r for b in blocks for r in [c[1] for c in b.cites] if r})
    single_only = sum(1 for b in blocks
                      if not any(audit.get(c[1]) == TWO_CARRIER for c in b.cites))
    with_proved = sum(1 for b in blocks
                      if any(ledger.get(c[1]) == PROVED for c in b.cites))
    print()
    print('  R11   SUMMARY  %d blocks · %d distinct rows cited · %d resting on no TWO-CARRIER row '
          '· %d carrying a PROVED row'
          % (len(blocks), len(distinct), single_only, with_proved))
    print()
    if fired:
        print('GATE FAILED — %d violation%s across %d rule%s'
              % (len(fired), '' if len(fired) == 1 else 's',
                 len(struck), '' if len(struck) == 1 else 's'))
        return 1
    print('GATE PASSED')
    return 0


if __name__ == '__main__':
    sys.exit(main())

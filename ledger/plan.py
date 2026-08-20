#!/usr/bin/env python3
"""THE PLAN — from where we are to a PROVEN process, retrieved and updated, never retyped.

  ./ledger/plan.py                 render the plan (and rewrite THE_PLAN_V001.md)
  ./ledger/plan.py set T-4 DONE    change one task's status
  ./ledger/plan.py add 4 T-20 "task text" TODO "done when..."   append a task

Same discipline as the status ledger: ledger/plan.tsv is the source of truth, the document is
generated, row order is file order, IDs never renumber, tasks append and are never deleted.
Every task carries a DONE_WHEN, so 'done' is checkable rather than a judgement call."""
import sys, os, csv

HERE = os.path.dirname(os.path.abspath(__file__)); ROOT = os.path.dirname(HERE)
PLAN = os.path.join(HERE, 'plan.tsv'); VOCAB = os.path.join(HERE, 'PLAN_VOCAB.tsv')
OUT = os.path.join(ROOT, 'THE_PLAN_V001.md')
MARK = {'DONE': '**DONE**', 'DOING': '**DOING**', 'TODO': 'TODO', 'BLOCKED': '**BLOCKED**', 'DROPPED': '~~DROPPED~~'}

def load(p):
    with open(p, newline='') as f: r = list(csv.reader(f, delimiter='\t'))
    return r[0], r[1:]
def esc(s): return s.replace('|', '\\|')

def render():
    hdr, rows = load(PLAN); _, vocab = load(VOCAB)
    valid = {v[0] for v in vocab}
    bad = [(r[0], r[3]) for r in rows if r[3] not in valid]
    if bad: sys.exit('status not in the closed vocabulary: %s' % bad)
    ids = [r[0] for r in rows]
    if len(ids) != len(set(ids)): sys.exit('duplicate task ids')

    done = sum(1 for r in rows if r[3] == 'DONE')
    L = ['# THE PLAN — to a PROVEN process', '',
         '**Call it with `./ledger/plan.py`. `ledger/plan.tsv` is the source of truth; this file is',
         'generated. Tasks append, IDs never renumber, and every task carries a DONE_WHEN so "done" is',
         'checkable rather than a judgement call.**', '',
         '## THE DEFINITION OF DONE', '',
         '> **A full demonstration of a PROVEN process that explains how EM, gravity, and alpha form',
         '> quantum records. And the model so that our proof can be checked by anyone.**',
         '>', '> — the principal, 2026-08-19', '',
         '**Nothing below counts as finished until both halves of that sentence are true.**', '',
         f'**{done} of {len(rows)} tasks done.**  ▶ marks the CRITICAL PATH — the chain the claim rests on.', '']
    tally = {}
    for r in rows: tally[r[3]] = tally.get(r[3], 0) + 1
    order = [s for s in ('DONE', 'DOING', 'TODO', 'BLOCKED', 'DROPPED') if s in tally]
    L += ['| ' + ' | '.join(order) + ' |', '|' + '---|' * len(order),
          '| ' + ' | '.join(str(tally[s]) for s in order) + ' |', '', '---', '']
    rows = sorted(rows, key=lambda r: int(r[7]) if len(r) > 7 and r[7] else 999)
    for ph in sorted({r[1] for r in rows}, key=lambda p: min(int(r[7]) for r in rows if r[1] == p)):
        sel = [r for r in rows if r[1] == ph]
        d = sum(1 for r in sel if r[3] == 'DONE')
        L += [f'## PHASE {ph}  —  {d}/{len(sel)} done', '',
              '| # | task | status | DONE WHEN | depends |', '|---|---|:---:|---|---|']
        for r in sel:
            i, _, task, st, dw, dep, row = (r + [''] * 7)[:7]
            crit = (len(r) > 8 and r[8] == 'yes')
            L.append('| %s**%s** | %s%s | %s | %s | %s |' %
                     ('▶ ' if crit else '', i, '**' if crit else '', esc(task) + ('**' if crit else ''),
                      MARK.get(st, st), esc(dw or '—'), esc(dep or '—')))
        L += ['']
    L += ['---', '', '## STATUS VOCABULARY', '', '| | |', '|---|---|']
    for v in vocab: L.append('| `%s` | %s |' % (v[0], v[1]))
    L += ['', '---', '',
          '> **The plan is not the goal. The program is. If a task turns out to be the wrong task,',
          '> it is DROPPED with the reason in the register, and a better one is appended.**', '']
    txt = '\n'.join(L)
    open(OUT, 'w').write(txt)
    return txt

def main():
    a = sys.argv[1:]
    if not a or a[0] == 'render': sys.stdout.write(render()); return
    hdr, rows = load(PLAN)
    if a[0] == 'set':
        _, vocab = load(VOCAB)
        if a[2] not in {v[0] for v in vocab}: sys.exit('not in the closed vocabulary: %s' % a[2])
        hit = [r for r in rows if r[0] == a[1]]
        if not hit: sys.exit('no such task: %s' % a[1])
        old = hit[0][3]; hit[0][3] = a[2]
    elif a[0] == 'add':
        # ORDER and CRIT were previously unreachable from the CLI, so appended tasks always sorted
        # last regardless of where they belong. add now accepts them: add PHASE ID TASK STATUS
        # [DONE_WHEN] [DEPENDS] [ROW] [ORDER] [CRIT]
        if any(r[0] == a[2] for r in rows): sys.exit('id already used')
        rows.append([a[2], a[1], a[3], a[4]] + (list(a[5:10]) + [''] * 5)[:5]); old = None
    elif a[0] == 'text':
        hit = [r for r in rows if r[0] == a[1]]
        if not hit: sys.exit('no such task: %s' % a[1])
        old = hit[0][2]; hit[0][2] = a[2]
    elif a[0] == 'dep':
        # DEPENDS is the plan's real logic -- a task blocked by another must be able to say so
        # without hand-editing the TSV, which the gate forbids.
        hit = [r for r in rows if r[0] == a[1]]
        if not hit: sys.exit('no such task: %s' % a[1])
        while len(hit[0]) < 9: hit[0].append('')
        old = hit[0][5]; hit[0][5] = a[2]
    else: sys.exit(__doc__)
    with open(PLAN, 'w', newline='') as f:
        w = csv.writer(f, delimiter='\t', lineterminator='\n'); w.writerow(hdr); w.writerows(rows)
    render(); print('%s: %s -> %s' % (a[1], old, a[2]) if old else 'appended %s' % a[2])

main()

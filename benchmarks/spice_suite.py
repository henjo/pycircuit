"""Run the public SPICE benchmark decks through pycircuit and score them
(the SPICE benchmark plan's stage 4).

    python benchmarks/spice_suite.py                 # every importable deck
    python benchmarks/spice_suite.py 4049osc gm1     # some
    python benchmarks/spice_suite.py --peers         # ngspice and Xyce too
    python benchmarks/spice_suite.py --json out.json
    python benchmarks/spice_suite.py --qpart 1 gm2   # the MOS gate charge's channel split

Each deck runs in a process of its own (its peak RSS is its own): read,
import, the `.tran` (`Imported.transient`), its wall times, steps, Newton
iterations and force-accepts, and each `.print tran` column scored by
Xyce's own metric (`spiceoutput.xyce_verify`, 1 or below passes) against
Xyce's gold where the suite has one.  `--peers` also writes the imported
circuit for ngspice (`Imported.write_ngspice`) and runs it, scoring
ngspice against the gold -- the three-way check that separates a model
difference (both away from the gold) from a simulator one -- or, where
there is no gold, scoring pycircuit against ngspice (an absolute floor of
1 % of each column's full scale, where a gold deck's columns carry
offsets); and times Xyce on the original deck.  The peers' files stay in the benchmark cache
(they derive from the decks, which declare no license); the suite never
calls a simulator.  Needs the fetched data (`benchmarks/fetch_spice_suite.py`).
"""
import json
import os
import re
import resource
import shutil
import subprocess
import sys
import tempfile
import time
import warnings

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from pycircuit._testing import benchdata

#: a deck's process is stopped after this many seconds, a peer's run after
#: `PEER_TIMEOUT` (their records say so; Xyce runs some decks for 10 min+)
CASE_TIMEOUT = 1800
PEER_TIMEOUT = 300
NGSPICE = os.environ.get('PYCIRCUIT_NGSPICE', os.path.expanduser('~/local/ngspice/bin/ngspice'))
XYCE = os.environ.get('PYCIRCUIT_XYCE', os.path.expanduser('~/local/xyce/serial/bin/Xyce'))

#: name -> (deck, gold or None): the decks with a `.tran` that import
#: (stage 2's census; latch and opampal since stage 7's substrate
#: junction) -- MOS level 2 (stage 9) keeps the rest out for now
_X = 'Netlists/'
_G = 'OutputData/'
CASES = {
    '4049osc': (_X + '4049OSC/4049osc.cir', _G + '4049OSC/4049osc.cir.prn'),
    'toronto': (_X + 'CircuitSim90/MOS2/toronto.cir', _G + 'CircuitSim90/MOS2/toronto.cir.prn'),
    'slowlatch': (_X + 'CircuitSim90/MOS2/slowlatch.cir',
                  _G + 'CircuitSim90/MOS2/slowlatch.cir.prn'),
    'rca': (_X + 'MCNC_BJT_RCA/rca.cir', _G + 'MCNC_BJT_RCA/rca.cir.prn'),
    'schmitecl': (_X + 'MCNC_BJT_SCHMITECL/schmitecl.cir_NORUN',
                  _G + 'MCNC_BJT_SCHMITECL/schmitecl.cir.prn'),
    'latch': (_X + 'MCNC_BJT_LATCH/latch.cir', _G + 'MCNC_BJT_LATCH/latch.cir.prn'),
    'opampal': (_X + 'MCNC_BJT_OPAMPAL/opampal.cir', _G + 'MCNC_BJT_OPAMPAL/opampal.cir.prn'),
}
for _n in ('gm1', 'gm2', 'gm3', 'gm17', 'mike2', 'rich3', 'todd3'):
    CASES[_n] = (f'{_X}CircuitSim90/MOS3/{_n}.cir', f'{_G}CircuitSim90/MOS3/{_n}.cir.prn')
for _n in ('arom', 'gm19', 'jge'):
    CASES[_n] = (f'{_X}CircuitSim90/MOS3/{_n}.cir', None)

#: a `.print` column: `v(a)`, `v(a,b)`, `i(vsrc)`, each perhaps `{... + c}`
_PROBE = re.compile(r'^\{?\s*([vi])\(\s*([^,()\s]+)\s*(?:,\s*([^()\s]+)\s*)?\)'
                    r'\s*(?:([+-])\s*([0-9.]+(?:[eE][-+]?\d+)?))?\s*\}?$', re.IGNORECASE)


def probes(net):
    """The `.print tran` columns: [(column, kind, nodes, offset)] -- those
    of another shape left out (named)."""
    out, skipped = [], []
    for d in net.prints:
        if not d.words or d.words[0].lower() != 'tran':
            continue
        for col in d.words[1:]:
            m = _PROBE.match(col)
            if m is None:
                skipped.append(col)
                continue
            kind, a, b, sign, c = m.groups()
            off = float(c) * (-1 if sign == '-' else 1) if c else 0.0
            out.append((col, kind.lower(), (a.lower(), b.lower() if b else None), off))
    return out, skipped


def _within(t, gold_t):
    """The test times inside the gold's span: a TSTOP read as SPICE reads it
    (`100n` is `100 * 1e-9`) can pass the gold's printed one by an ulp --
    within 1e-12 of the span it is clamped, beyond it dropped."""
    import numpy as np
    span = gold_t[-1] - gold_t[0]
    keep = (t >= gold_t[0] - 1e-12 * span) & (t <= gold_t[-1] + 1e-12 * span)
    return np.clip(t[keep], gold_t[0], gold_t[-1]), keep


def _gold_column(names, col):
    lower = [n.lower() for n in names]
    key = col.lower()
    return lower.index(key) if key in lower else None


def run_one(name, peers, qpart=None):
    """One deck, in this process: its record (see the module note)."""
    import numpy as np

    from pycircuit.circuit import circuit, spice_import
    from pycircuit.utilities import spicenetlist, spiceoutput
    deck_rel, gold_rel = CASES[name]
    deck = benchdata.spice_data(deck_rel)
    gold = benchdata.spice_data(gold_rel) if gold_rel else None
    if deck is None or (gold_rel and gold is None):
        return {'case': name, 'skipped': 'data not fetched'}
    circuit.default_toolkit = circuit.numeric
    rec = {'case': name}
    t0 = time.perf_counter()
    net = spicenetlist.read(deck)
    rec['read_s'] = time.perf_counter() - t0
    t0 = time.perf_counter()
    imp = spice_import.import_netlist(deck, mos_qpart=qpart)
    rec['import_s'] = time.perf_counter() - t0
    rec['elements'], rec['unknowns'] = len(imp.elements), imp.circuit.n
    tr, kw = imp.transient()
    t0 = time.perf_counter()
    with warnings.catch_warnings():
        warnings.simplefilter('ignore')
        res = tr.solve(**kw)
    rec['solve_s'] = time.perf_counter() - t0
    st = tr.statistics
    rec.update(steps=st.accepted_steps, rejected=st.rejected_steps,
               newton=st.newton_iterations, force_accepts=st.force_accepts)
    rec['rss_mb'] = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 1024
    cols, skipped = probes(net)
    rec['skipped_columns'] = skipped
    comps = spiceoutput.comp_tolerances([c for c, _w in net.comps])
    t = None
    ours = {}
    for col, kind, (a, b), off in cols:
        w = res.v(a, b) if kind == 'v' and b else res.v(a) if kind == 'v' else res.i(f'{a}.plus')
        t = np.asarray(w.x[0], dtype=float)
        ours[col] = np.asarray(w.y, dtype=float) + off
    rec['metric'] = {}
    if gold is not None:
        names, g = spiceoutput.read_prn(gold)
        ts, keep = _within(t, g[:, 1])
        for col in ours:
            k = _gold_column(names, col)
            if k is not None:
                rec['metric'][col] = spiceoutput.xyce_verify(ts, ours[col][keep], g[:, 1],
                                                             g[:, k], **comps.get(col, {}))
    if peers:
        out, ng = _peers(name, imp, deck, cols, gold, comps)
        rec.update(out)
        if gold is None and ng is not None:
            ## (no gold: ngspice's waveform is the reference -- with an
            ## absolute floor of 1 % of its full scale, where a gold deck's
            ## columns carry offsets instead: a logic low at 0 V would
            ## otherwise read any difference as relative error)
            nt, nv = ng
            ts, keep = _within(t, nt)
            rec['metric_vs_ngspice'] = {
                col: spiceoutput.xyce_verify(
                    ts, ours[col][keep], nt, nv[col],
                    **{'abstol': 1e-2 * float(np.max(np.abs(nv[col]))), **comps.get(col, {})})
                for col in ours if col in nv}
    return rec


def _peers(name, imp, deck, cols, gold, comps):
    """ngspice on the imported circuit written back, scored against the
    gold; Xyce's wall time on the original deck: `(record, (ngspice's
    times, {column: its waveform}) or None)`."""
    from pycircuit.utilities import spiceoutput
    out, ng = {}, None
    work = os.path.join(benchdata.cache_root(), 'runs', name)
    os.makedirs(work, exist_ok=True)
    nodes = sorted({n for _c, kind, pair, _o in cols if kind == 'v' for n in pair if n})
    currents = sorted({('i', pair[0]) for _c, kind, pair, _o in cols if kind == 'i'})
    path = os.path.join(work, f'{name}.ng.cir')
    names = imp.write_ngspice(path, probes=nodes + currents)
    if shutil.which(NGSPICE) or os.path.exists(NGSPICE):
        t0 = time.perf_counter()
        try:
            p = subprocess.run([NGSPICE, '-b', path], capture_output=True, text=True,
                               cwd=work, timeout=PEER_TIMEOUT, check=False)
        except subprocess.TimeoutExpired:
            out['ngspice_s'] = time.perf_counter() - t0
            out['ngspice_error'] = f'over {PEER_TIMEOUT} s'
            return out, None
        out['ngspice_s'] = time.perf_counter() - t0
        try:
            ncols, data = spiceoutput.read_ngspice_print(p.stdout)
        except ValueError as e:
            out['ngspice_error'] = str(e)
        else:
            low = [c.lower() for c in ncols]
            out['ngspice_steps'] = len(data) - 1
            waves = {}
            for col, kind, (a, b), off in cols:
                if kind == 'i':
                    va, vb = data[:, low.index(names[('i', a)].lower())], 0.0
                else:
                    va = data[:, low.index(f'v({names[a]})')]
                    vb = data[:, low.index(f'v({names[b]})')] if b else 0.0
                waves[col] = va - vb + off
            ng = (data[:, 0], waves)
            if gold is not None:
                gn, g = spiceoutput.read_prn(gold)
                out['ngspice_metric'] = {}
                ts, keep = _within(data[:, 0], g[:, 1])
                for col, wave in waves.items():
                    k = _gold_column(gn, col)
                    if k is not None:
                        out['ngspice_metric'][col] = spiceoutput.xyce_verify(
                            ts, wave[keep], g[:, 1], g[:, k], **comps.get(col, {}))
    if shutil.which(XYCE) or os.path.exists(XYCE):
        with tempfile.TemporaryDirectory(dir=work) as tmp:
            local = os.path.join(tmp, os.path.basename(deck))
            shutil.copy(deck, local)
            t0 = time.perf_counter()
            try:
                p = subprocess.run([XYCE, local], capture_output=True, text=True, cwd=tmp,
                                   timeout=PEER_TIMEOUT, check=False)
                out['xyce_ok'] = p.returncode == 0
            except subprocess.TimeoutExpired:
                out['xyce_ok'] = False
            out['xyce_s'] = time.perf_counter() - t0
    return out, ng


def _fmt(rec):
    if 'skipped' in rec:
        return f'{rec["case"]:10s} skipped: {rec["skipped"]}'
    worst = max(rec['metric'].values()) if rec['metric'] else None
    s = (f'{rec["case"]:10s} {rec["elements"]:6d} el {rec["unknowns"]:6d} unk  '
         f'read {rec["read_s"]:6.2f}s import {rec["import_s"]:6.2f}s solve {rec["solve_s"]:7.2f}s  '
         f'{rec["steps"]:6d} steps {rec["rejected"]:5d} rej {rec["newton"]:7d} newton '
         f'{rec["force_accepts"]:3d} forced  {rec["rss_mb"]:6.0f} MB  '
         + ('metric ' + (f'{worst:.3f}' if worst is not None else '  -  ')))
    if rec.get('metric_vs_ngspice'):
        s += f' (vs ngspice: {max(rec["metric_vs_ngspice"].values()):.3f})'
    if 'ngspice_s' in rec:
        nm = rec.get('ngspice_metric') or {}
        s += (f'  | ngspice {rec["ngspice_s"]:6.2f}s'
              + (f' metric {max(nm.values()):.3f}' if nm else ''))
    if 'xyce_s' in rec:
        s += f'  | Xyce {rec["xyce_s"]:6.2f}s' + ('' if rec.get('xyce_ok') else ' (failed)')
    return s


def main(argv):
    peers = '--peers' in argv
    qpart = float(argv[argv.index('--qpart') + 1]) if '--qpart' in argv else None
    out_json = argv[argv.index('--json') + 1] if '--json' in argv else None
    if '--one' in argv:
        print(json.dumps(run_one(argv[argv.index('--one') + 1], peers, qpart)))
        return 0
    values = {out_json, argv[argv.index('--qpart') + 1] if '--qpart' in argv else None}
    names = [a for a in argv if not a.startswith('--') and a not in values] or list(CASES)
    unknown = [n for n in names if n not in CASES]
    if unknown:
        sys.exit(f'unknown cases {unknown}; known: {", ".join(CASES)}')
    recs = []
    for name in names:
        cmd = [sys.executable, os.path.abspath(__file__), '--one', name] + (
            ['--peers'] if peers else []) + (['--qpart', str(qpart)] if qpart is not None else [])
        try:
            p = subprocess.run(cmd, capture_output=True, text=True, check=False,
                               timeout=CASE_TIMEOUT)
        except subprocess.TimeoutExpired:
            rec = {'case': name, 'skipped': f'over {CASE_TIMEOUT} s'}
            recs.append(rec)
            print(_fmt(rec), flush=True)
            continue
        line = p.stdout.strip().splitlines()[-1] if p.stdout.strip() else ''
        try:
            rec = json.loads(line)
        except ValueError:
            rec = {'case': name, 'skipped': 'failed: ' + (p.stderr.strip().splitlines() or
                                                          ['?'])[-1][:200]}
        recs.append(rec)
        print(_fmt(rec), flush=True)
    if out_json:
        with open(out_json, 'w') as fh:
            json.dump(recs, fh, indent=1)
    return 0


if __name__ == '__main__':
    sys.exit(main(sys.argv[1:]))

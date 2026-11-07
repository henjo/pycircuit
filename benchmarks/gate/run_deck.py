"""One gold deck, one method, one setting: `run_deck.py deck method setting` -> one JSON line.  Run from the tree under test (cwd and
PYTHONPATH): `file` says which package it imported.  The paths on or off are the environment's (`settings.PATHS`)."""
import hashlib
import json
import logging
import sys
import time
import warnings

warnings.simplefilter('ignore')
logging.disable(logging.WARNING)
sys.path.insert(0, 'benchmarks')
import numpy as np
import settings
import spice_suite

import pycircuit
from pycircuit.circuit import _paths, circuit, spice_import
from pycircuit.circuit import integrator as I

KLU = ('factors', 'refactors', 'pivot_fallbacks', 'residual_fallbacks', 'chord_factors', 'chord_refactors')


def _sha(*arrays):
    h = hashlib.sha1()
    for a in arrays:
        h.update(np.asarray(a, float).tobytes())
    return h.hexdigest()[:12]


def _swing(x, t):
    """How many unknowns move over 100 mV: a benchmark that barely moves times and tests little."""
    x = np.asarray(x, float)
    s = x.max(axis=0 if x.shape[0] == np.size(t) else -1) - x.min(axis=0 if x.shape[0] == np.size(t) else -1)
    return int((s > 0.1).sum())


def _dc(imp, rec):
    from pycircuit.circuit.dcanalysis import DC
    tr, _kw = imp.transient()
    ns = dict(tr.par.nodeset) if getattr(tr.par, 'nodeset', None) else None
    dc = DC(imp.circuit, epar=imp.epar(), nodeset=ns)
    first = dc.solve()
    second = DC(imp.circuit, epar=imp.epar(), nodeset=ns).solve()
    rec.update(x=_sha(first.x, second.x))
    return dc._get_linearsolver()


def _tran(imp, deck, method, setting, rec):
    tr, kw = imp.transient(**spice_suite.OVERRIDES.get(deck, {}))
    if method == 'radau':
        tr.par.integrator = I.RadauIIA3Integrator()
    settings.force_run(setting, tr)
    try:
        res = tr.solve(**kw)
        rec.update(x=_sha(res.x, res.sweep_values), swing=_swing(res.x, res.sweep_values))
    except Exception as e:  # noqa: BLE001 -- a run's failure is a result both paths must share
        rec.update(x='raised', err=type(e).__name__ + ': ' + str(e)[:100])
    st = tr.statistics
    rec.update(steps=st.accepted_steps, rej=st.rejected_steps, newton=st.newton_iterations, fb=getattr(st, 'chord_fallbacks', 0))
    return tr._get_linearsolver()


def main():
    deck, method, setting = sys.argv[1:4]
    circuit.default_toolkit = circuit.numeric
    settings.force_modules(setting)
    imp = spice_import.import_netlist(spice_suite.benchdata.spice_data(spice_suite.CASES[deck][0]))
    rec = {'deck': deck, 'method': method, 'setting': setting, 'n': imp.circuit.n, 'file': pycircuit.__file__}
    b = _paths.snapshot()
    t0 = time.perf_counter()
    ls = _dc(imp, rec) if method == 'dc' else _tran(imp, deck, method, setting, rec)
    rec['s'] = round(time.perf_counter() - t0, 3)
    ch = getattr(ls, '_choice', ls)
    rec.update(solver=type(ch).__name__, klu={k: getattr(ch, k, None) for k in KLU}, counts=_paths.since(b))
    print(json.dumps(rec))


if __name__ == '__main__':
    main()

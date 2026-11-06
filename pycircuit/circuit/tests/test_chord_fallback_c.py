"""The chord's fallback in the C Newton (newton.h's THE CHORD'S FALLBACK,
`_tran_newton_c.FALLBACK`; roadmap 268).  `ChordNewton.solve_system` falls
back to the full Newton from the same seed where its chord breaks or meets
its limit; the C Newton handed those attempts back, and the Python repeated
the whole attempt on dense state -- every Gear-2 hand-back measured was one
(a 4 000-stage chain's 1 us run: 12, about 70 % of its time).  Now the C
runs the full Newton there: its first evaluation the chord's first C (what
`_C_lookup` finds), no early check past the chord's, its iterations added,
`chord_fallbacks` counted, a sparse run's numeric refactored as the
Python's full Newton refactors it -- every result the Python's."""
import warnings

import numpy as np
import pytest

from pycircuit._testing import benchdata
from pycircuit.circuit import (
    _paths,
    _tran_loop_c,
    _tran_lte_c,
    _tran_newton_c,
    circuit,
    spice_import,
)

PATHS = {'python': (False, False, False), 'python loop, C Newton': (False, True, True),
         'C loop': (True, True, True)}


def _deck(monkeypatch, name, path, fallback=True, dense=False):
    """todd3 / slowlatch (CircuitSim90, the suite's way) on `path`: the run's
    bytes and statistics, its linear solver's refactors, its counts."""
    data = benchdata.spice_data(f'Netlists/CircuitSim90/{name}.cir')
    if data is None:
        pytest.skip('benchmark data not fetched (benchmarks/fetch_spice_suite.py)')
    loop, newton, lte = PATHS[path]
    monkeypatch.setattr(_tran_loop_c, 'ENABLED', loop)
    monkeypatch.setattr(_tran_newton_c, 'ENABLED', newton)
    monkeypatch.setattr(_tran_lte_c, 'ENABLED', lte)
    monkeypatch.setattr(_tran_newton_c, 'FALLBACK', fallback)
    circuit.default_toolkit = circuit.numeric
    imp = spice_import.import_netlist(data)
    tr, kw = imp.transient()
    if dense:
        ## (numpy's dense solve whatever the size: the full Newton's dgesv)
        ls = tr._get_linearsolver()
        ls.min_n = ls.structural_min_n = 10 ** 9
    before = _paths.snapshot()
    with warnings.catch_warnings():
        warnings.simplefilter('ignore')
        res = tr.solve(**kw)
    st = tr.statistics
    ch = tr._get_linearsolver()._choice
    run = (np.asarray(res.x, float).tobytes(), np.asarray(res.sweep_values, float).tobytes(),
           st.accepted_steps, st.rejected_steps, st.newton_iterations, st.chord_fallbacks,
           st.early_aborts)
    klu = (type(ch).__name__, getattr(ch, 'refactors', None), getattr(ch, 'factors', None))
    return run, klu, _paths.since(before)


@pytest.mark.parametrize('dense', [False, True], ids=['sparse', 'dense'])
@pytest.mark.parametrize('name', ['MOS3/todd3', 'MOS2/slowlatch'])
def test_the_c_falls_back_as_the_python_does(monkeypatch, name, dense):
    """todd3 (chord breaks and limits) and slowlatch (a limit), on KLU and on
    numpy's dense solve: the C loop's run with the fallback is the Python's
    -- bytes, steps, rejections, Newton iterations, chord fallbacks, early
    aborts, the KLU's refactors and factors -- and the C loop's without it
    (the attempts handed back); no chord hands an attempt back."""
    want, wklu, wcounts = _deck(monkeypatch, name, 'python', dense=dense)
    got, gklu, counts = _deck(monkeypatch, name, 'C loop', dense=dense)
    off, oklu, ocounts = _deck(monkeypatch, name, 'C loop', fallback=False, dense=dense)
    assert want[5] > 0 and wklu[0] == ('DenseSolver' if dense else 'StructuralKLU')
    assert got == want == off
    assert gklu == wklu == oklu
    assert counts['newton_c:fallbacks'] == want[5]
    assert counts.get('early:checks', 0) == wcounts.get('early:checks', 0) > 0
    assert counts.get('loop_c:bail:newton', 0) == 0
    assert not any(k.startswith('newton_c:bail:') for k in counts)
    assert ocounts['loop_c:bail:newton'] > 0 and 'newton_c:fallbacks' not in ocounts


def test_the_python_loops_c_newton_falls_back_too(monkeypatch):
    """The Python loop on the C Newton (`_tran_newton_c.solve`): the same
    run, the fallbacks served, none handed back."""
    want, wklu, wcounts = _deck(monkeypatch, 'MOS3/todd3', 'python')
    got, gklu, counts = _deck(monkeypatch, 'MOS3/todd3', 'python loop, C Newton')
    assert got == want and gklu == wklu
    assert counts.get('early:checks', 0) == wcounts.get('early:checks', 0) > 0
    assert counts['newton_c:fallbacks'] == want[5] > 0
    assert not any(k.startswith('newton_c:bail:') for k in counts)

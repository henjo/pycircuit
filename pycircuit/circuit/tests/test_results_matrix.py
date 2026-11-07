"""A RUN'S RESULTS HELD ONCE (`_tran_history.results_matrix`; plan 261's L14):
the run's accepted states, one array a step, become one matrix.  numpy's copy
holds both at its end; a large run's rows are copied one at a time, each
released and the C heap trimmed as the matrix grows.  The matrix is numpy's to
the bit, and what cannot be gathered that way is left to numpy."""
import subprocess
import sys
import textwrap
import types

import numpy as np
import pytest

from pycircuit.circuit import _tran_history as th
from pycircuit.circuit import circuit
from pycircuit.circuit.integrator import RadauIIA3Integrator
from pycircuit.circuit.tests.test_hdl_batch import mos_chain
from pycircuit.circuit.transient import Transient


def _rows():
    """Rows holding what a copy could disturb: a NaN's payload and sign, the
    signed zeros, the subnormals, the infinities."""
    rng = np.random.default_rng(5)
    X = [rng.standard_normal(9) for _ in range(40)]
    nan = np.frombuffer(np.array([0xFFF8_0000_DEAD_BEEF], dtype=np.uint64).tobytes(), dtype=np.float64)[0]
    X[3][2] = nan
    X[7][:4] = -0.0, 5e-324, np.inf, -np.inf
    return X


def test_small_results_are_numpys_own_and_the_list_stays(monkeypatch):
    """Under the size nothing is gathered: the list as it was, numpy's array."""
    X = _rows()
    kept = list(X)
    monkeypatch.setattr(th, '_malloc_trim', lambda: pytest.fail('nothing to trim under the size'))
    M = th.results_matrix(X)
    assert len(X) == len(kept) and all(a is b for a, b in zip(X, kept))
    assert M.tobytes() == np.array(kept).tobytes()


def test_large_results_are_the_same_bytes_and_the_rows_are_released(monkeypatch):
    """From the size the matrix is numpy's to the bit -- a NaN's payload, the
    signed zeros, the subnormals -- its shape, type and layout too, and the
    list is spent."""
    X = _rows()
    want = np.array(X)
    monkeypatch.setattr(th, 'GATHER_MIN_BYTES', 1)
    calls = []
    monkeypatch.setattr(th, '_malloc_trim', lambda: calls.append)
    M = th.results_matrix(X)
    assert X == []
    assert M.tobytes() == want.tobytes()
    assert M.shape == want.shape and M.dtype == want.dtype and M.flags.c_contiguous and M.flags.writeable
    assert calls, 'the heap is trimmed'


def test_no_malloc_trim_gathers_all_the_same(monkeypatch):
    """A C library without `malloc_trim` gathers without trimming."""
    X = _rows()
    want = np.array(X).tobytes()
    monkeypatch.setattr(th, 'GATHER_MIN_BYTES', 1)
    monkeypatch.setattr(th, '_malloc_trim', lambda: None)
    assert th.results_matrix(X).tobytes() == want and X == []


@pytest.mark.parametrize('rows', [
    [np.arange(3.0), np.arange(3.0) * 1j],
    [np.arange(3), np.arange(3)],
    [np.arange(3.0), [1.0, 2.0, 3.0]],
    [np.arange(3.0).reshape(3, 1), np.arange(3.0).reshape(3, 1)],
], ids=['complex', 'integer', 'a-list', 'two-dimensional'])
def test_rows_that_are_not_float64_arrays_are_numpys(monkeypatch, rows):
    """Another type, a row that is not an array, another shape: numpy's array,
    the list as it was."""
    monkeypatch.setattr(th, 'GATHER_MIN_BYTES', 1)
    X = list(rows)
    M = th.results_matrix(X)
    assert len(X) == len(rows)
    want = np.array(rows)
    assert M.dtype == want.dtype and M.tobytes() == want.tobytes()


def test_rows_of_unequal_length_are_refused_as_numpy_refuses(monkeypatch):
    monkeypatch.setattr(th, 'GATHER_MIN_BYTES', 1)
    X = [np.arange(3.0), np.arange(4.0)]
    with pytest.raises(ValueError):
        np.array(X)
    with pytest.raises(ValueError):
        th.results_matrix(X)
    assert len(X) == 2


def test_no_states_are_numpys(monkeypatch):
    monkeypatch.setattr(th, 'GATHER_MIN_BYTES', 1)
    assert th.results_matrix([]).shape == (0,)


def test_a_large_gathering_holds_the_rows_once():
    """A fresh process's peak over the rows' own: under half the results
    (numpy's copy, or releasing without trimming: all of them)."""
    pytest.importorskip('ctypes')
    if th._malloc_trim() is None:
        pytest.skip('no malloc_trim here')
    code = textwrap.dedent('''
        import numpy as np
        from pycircuit.circuit import _tran_history as th
        def hwm():
            for l in open('/proc/self/status'):
                if l.startswith('VmHWM'):
                    return int(l.split()[1]) * 1024
        n, rows = 4005, 10000
        Xo = np.random.default_rng(1).random(n)
        X = [Xo.copy() for _ in range(rows)]
        before = hwm()
        M = th.results_matrix(X)
        print((hwm() - before) / M.nbytes, M.shape == (rows, n), float(M[-1, -1]) == float(Xo[-1]))
    ''')
    out = subprocess.run([sys.executable, '-c', code], capture_output=True, text=True, timeout=120, check=False)
    assert out.returncode == 0, out.stderr
    rise, shape, last = out.stdout.split()
    assert shape == 'True' and last == 'True'
    assert float(rise) < 0.5, rise


@pytest.mark.parametrize('radau', [False, True], ids=['gear2', 'radau'])
def test_a_run_s_result_is_the_same_gathered(monkeypatch, radau):
    """A whole run: gathered or not its result's bytes are the same, its list
    spent, what it was made from the run's own."""
    def run():
        kw = {'integrator': RadauIIA3Integrator()} if radau else {}
        tr = Transient(mos_chain(6), toolkit=circuit.numeric, **kw)
        return tr, tr.solve(tend=2e-8, timestep=1e-9)
    _tr, plain = run()
    seen = []
    real = th.results_matrix

    def spy(X):
        seen.append(X)
        return real(X)
    monkeypatch.setattr(th, 'GATHER_MIN_BYTES', 1)
    from pycircuit.circuit import transient as T
    monkeypatch.setattr(T, 'results_matrix', spy)
    tr, got = run()
    assert len(seen) == 1 and seen[0] == []
    assert np.asarray(got.x).tobytes() == np.asarray(plain.x).tobytes()
    assert np.asarray(got.sweep_values).tobytes() == np.asarray(plain.sweep_values).tobytes()
    assert got.statistics is tr.statistics


def test_another_toolkit_makes_its_own_array(monkeypatch):
    """A toolkit other than the numeric's makes the matrix with its `array`."""
    tr = Transient(mos_chain(4), toolkit=circuit.numeric)
    tr.solve(tend=3e-9, timestep=1e-9)
    from pycircuit.circuit import transient as T
    monkeypatch.setattr(T, 'results_matrix', lambda X: pytest.fail('the numeric toolkit\'s gathering'))
    made = []

    def array(x):
        made.append(x)
        return np.array(x)
    tr.toolkit = types.SimpleNamespace(array=array)
    X = [np.arange(3.0), np.arange(3.0) + 1.0]
    res = tr._finish_result(X, [1e-9], 0.0)
    assert made[0] is X and np.asarray(res.x).shape == (3, 2)

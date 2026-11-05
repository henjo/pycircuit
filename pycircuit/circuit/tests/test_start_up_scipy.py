"""A run's start-up without scipy: the package imports none of it, and a
transient or DC of the gold decks reaches none on its served path.

scipy.sparse (+ .linalg) was 0.12 s of every process's import and
scipy.linalg 0.10 s more of most runs (2026-11-04): the sparse toolkit's
backend and the PSS walks defer it (`_lazy.deferred`), the branch check's
connected components are numpy's (`_tran_branch._components`), the C Newton
maps scipy's OpenBLAS by its path (`_cbuild.openblas('scipy')`), and the
complex factor's product loads scipy's own `_sparsetools` alone
(`linearsolver._scipy_sparsetools`) -- each the same code or the same answer.
"""
import ctypes
import glob
import os
import subprocess
import sys
import textwrap

import numpy as np

ROOT = os.path.dirname(os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__)))))


def _python(code):
    out = subprocess.run([sys.executable, '-c', textwrap.dedent(code)], cwd=ROOT,
                         capture_output=True, text=True, timeout=900, check=False)
    assert out.returncode == 0, out.stderr[-4000:]
    return out.stdout.strip().splitlines()


_DECKS = '''
    import sys, warnings, logging
    warnings.simplefilter('ignore'); logging.disable(logging.WARNING)
    sys.path.insert(0, 'benchmarks')
    import spice_suite
    from pycircuit.circuit import spice_import, circuit, integrator
    from pycircuit.circuit.dcanalysis import DC
    from pycircuit.circuit.transient import Transient
    assert not any(m == 'scipy' or m.startswith('scipy.') for m in sys.modules)
    circuit.default_toolkit = circuit.numeric
    for deck in ('rca', 'gm19', 'voter25'):
        for what in ('gear2', 'radau', 'dc'):
            imp = spice_import.import_netlist(
                spice_suite.benchdata.spice_data(spice_suite.CASES[deck][0]))
            tr, kw = imp.transient(**spice_suite.OVERRIDES.get(deck, {}))
            if what == 'dc':
                ns = dict(tr.par.nodeset) if getattr(tr.par, 'nodeset', None) else None
                DC(imp.circuit, epar=imp.epar(), nodeset=ns).solve()
                continue
            if what == 'radau':
                tr.par.integrator = integrator.RadauIIA3Integrator()
            tr.solve(**kw)
    print(sorted(m for m in sys.modules if m == 'scipy' or m.startswith('scipy.')))
'''


def test_a_warm_run_of_three_decks_imports_no_scipy():
    """Three gold decks -- a Gummel-Poon bipolar and two large MOS -- under
    gear2, Radau and DC by a second process (the first may compile): no
    scipy module is imported."""
    _python(_DECKS)
    assert _python(_DECKS)[-1] == '[]'


def test_the_components_are_scipy_s():
    """`_components` labels as `connected_components(directed=True,
    connection='weak')` does -- count, labels, their dtype -- on random
    patterns, chains and the empty graph."""
    from scipy.sparse import csr_matrix
    from scipy.sparse.csgraph import connected_components

    from pycircuit.circuit._tran_branch import _components
    rng = np.random.default_rng(7)
    cases = [np.zeros((1, 1), bool), np.zeros((5, 5), bool)]
    for _ in range(400):
        n = int(rng.integers(1, 50))
        cases.append(rng.random((n, n)) < rng.choice([0.01, 0.04, 0.1, 0.4]))
        chain = np.zeros((n, n), bool)
        idx = rng.permutation(n)
        chain[idx[:-1], idx[1:]] = rng.random(n - 1) < 0.85
        cases.append(chain)
    for pat in cases:
        n = pat.shape[0]
        r, c = np.nonzero(pat)
        want = connected_components(csr_matrix((np.ones(r.size, np.int8), (r, c)), shape=(n, n)),
                                    directed=True, connection='weak')
        got = _components(n, r, c)
        assert got[0] == want[0] and got[1].dtype == want[1].dtype
        assert got[1].tobytes() == want[1].tobytes()


def test_the_product_loaded_alone_is_scipy_s():
    """The complex factor's product: `_sparsetools` loaded alone gives
    scipy's `csc_matrix.dot` bytes, and where scipy.sparse is imported it is
    scipy's own module."""
    import scipy.sparse as S

    from pycircuit.circuit import linearsolver as ls
    alone = None
    saved = sys.modules.pop('scipy.sparse._sparsetools')
    try:
        alone = ls._scipy_sparsetools()
    finally:
        sys.modules['scipy.sparse._sparsetools'] = saved
    assert alone is not saved and alone.__name__ == 'pycircuit._scipy._sparsetools'
    assert ls._scipy_sparsetools() is saved
    rng = np.random.default_rng(3)
    for _ in range(100):
        n = int(rng.integers(1, 60))
        A = (rng.random((n, n)) < 0.15) * (rng.standard_normal((n, n)) + 1j * rng.standard_normal((n, n)))
        m = S.csc_matrix(A).astype(np.complex128)
        x = rng.standard_normal(n) + 1j * rng.standard_normal(n)
        got = ls._csc_dot(alone.csc_matvec, n, m.indptr, m.indices, m.data, x)
        assert got.tobytes() == m.dot(x).tobytes()


def test_scipy_s_openblas_by_name_is_the_one_scipy_maps():
    """`openblas('scipy')` (scipy not imported for it) is the library scipy's
    LAPACK runs: one mapping of the file, the routines at the addresses a
    fresh load of its path gives."""
    import scipy
    import scipy.linalg

    from pycircuit.circuit import _cbuild
    lib = _cbuild.openblas('scipy')
    assert lib is _cbuild.openblas(scipy)
    with open('/proc/self/maps') as fh:
        maps = {ln.split()[-1] for ln in fh if 'libscipy_openblas-' in ln}
    assert len(maps) == 1
    d = os.path.join(os.path.dirname(scipy.__file__), '..', 'scipy.libs')
    again = ctypes.CDLL(os.path.realpath(glob.glob(os.path.join(d, 'libscipy_openblas-*.so'))[0]))
    assert _cbuild.address(lib.scipy_dgetrf_) == _cbuild.address(again.scipy_dgetrf_)


def test_the_sparse_toolkit_imports_its_backend_on_first_use():
    """The sparse toolkit's backend is deferred: a fresh process imports
    the package without scipy, and the toolkit's first primitive brings it
    in and works."""
    lines = _python('''
        import sys
        from pycircuit.circuit import toolkit
        assert 'scipy' not in sys.modules and 'pycircuit.circuit._sparse_numeric' not in sys.modules
        z = toolkit.sparse_numeric.zeros((2, 2))
        print(type(z).__name__, 'pycircuit.circuit._sparse_numeric' in sys.modules)
    ''')
    assert lines[-1].endswith(' True')


def test_the_package_maps_scipy_s_openblas_without_scipy():
    """The package's import maps scipy's OpenBLAS (its threads start over the
    start-up, not mid-solve, where their spin cost every run ~0.1 s) and
    imports no scipy module."""
    lines = _python('''
        import sys
        import pycircuit.circuit
        with open('/proc/self/maps') as fh:
            mapped = any('libscipy_openblas-' in ln for ln in fh)
        print(mapped, any(m == 'scipy' or m.startswith('scipy.') for m in sys.modules))
    ''')
    assert lines[-1] == 'True False'

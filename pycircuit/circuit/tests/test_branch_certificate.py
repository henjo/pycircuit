"""THE BRANCH CHECK'S CERTIFICATE (`_BranchCheck._branch_certify`,
2026-10-14; Andreas: "Build it properly"; doc/hdl_roadmap_260824.md, 133).

Where the screen's cheap level fires, `rank C >= r0` is proved before the
SVD: sigma_r0(C) >= sigma_min(B) for an r0 x r0 submatrix B
(interlacing), B block diagonal, sigma_min(B_k) >= 1 / ||B_k^-1||_F.
Measured on the 27 gold decks: four fired the cheap level on every step
and the SVD never fired (0 of 5 834) -- most of their time, and it kept them
out of the C loop; the certificate proved all 5 834.

What these hold: a proof is the SVD's verdict without the SVD; what it
cannot prove the SVD decides -- a rank that holds with a bound below tol,
an entry between two blocks, a submatrix that is singular; a rank that
drops is never proved and still fires; the C's certificate is the Python's
bit for bit, across block sizes about OpenBLAS's thresholds and margins
about tol, threaded and on one thread; a run is the same with it on and
off on every path, and the C loop takes the attempts the screen handed
back; its data are made once, with the structural rank."""
import warnings

import numpy as np
import pytest
from threadpoolctl import threadpool_limits

from pycircuit._testing import benchdata
from pycircuit.circuit import _paths, _tran_branch, _tran_branch_c, _tran_loop_c, circuit
from pycircuit.circuit._tran_branch import _BranchCheck, _certify
from pycircuit.circuit.elements import C, R, SubCircuit, VPulse, gnd
from pycircuit.circuit.transient import Transient

pytestmark = pytest.mark.skipif(_tran_branch_c.driver() is None,
                                reason='the branch screen unit did not build')


class _Check(_BranchCheck):
    """The screen and its certificate alone: the run's structural rank and
    generic C set by hand (`_branch_structural_rank`'s), the SVD counted."""

    def __init__(self, C0, r0, iref=0):
        self._branch_rank0 = (r0, 1.0)
        self._branch_generic = (C0, C0 != 0)
        self.irefnode = iref
        self.svds = 0

    def _branch_svd(self, C, r0, ref):
        self.svds += 1
        return super()._branch_svd(C, r0, ref)


def _python(C, C0, r0):
    """`_branch_screen_C`'s verdict, and the SVDs it made."""
    t = _Check(C0, r0)
    verdict, _nrm = t._branch_screen_C(C, r0)
    return verdict, t.svds


def _c(C, C0, r0, pat=None):
    """The C screen's code (QUIET, COLLAPSED, SVD), its certificate's
    smallest bound where it proved, and its `branch_c:cert` count."""
    t = _Check(C0, r0)
    if pat is not None:
        t._C_pattern = (C, pat)
    before = _paths.snapshot()
    got = _tran_branch_c.screen(t, C)
    assert got is not None, _paths.since(before)
    drv = _tran_branch_c.driver()
    return got[0], (drv[2].bound if got[0] == _tran_branch_c.QUIET else None), \
        _paths.since(before).get('branch_c:cert', 0)


def _blocks(B3, big=1.0):
    """C (5 x 5): row 0 the reference, unknown 1 a 1 x 1 block `big` (the
    running maximum: tol = 1e-9 big), unknowns 2..4 the 3 x 3 block `B3`."""
    Cm = np.zeros((5, 5))
    Cm[1, 1] = big
    Cm[2:, 2:] = B3
    return Cm


#: an orthogonal 3 x 3 (fixed): s Q has every singular value s
_Q = np.linalg.qr(np.array([[2.0, 1.0, 0.5], [0.3, 1.5, 1.0], [1.0, 0.2, 1.7]]))[0]


def _sv(s):
    """A 3 x 3 with singular values `s`."""
    U = np.linalg.qr(np.array([[1.0, 0.4, 0.2], [0.1, 1.0, 0.7], [0.5, 0.3, 1.0]]))[0]
    return U @ np.diag(s) @ _Q.T


def test_a_proof_is_the_svds_verdict_without_the_svd():
    """Every singular value 1e-6 (tol 1e-9; the cheap level fires on the
    diagonal under 1e-6 of the maximum): quiet, no SVD, and the C proves
    it too -- on the 3 x 3 block's dgesv."""
    Cm = _blocks(1e-6 * _Q)
    assert np.abs(np.diag(Cm)[2:]).min() < 1e-6
    assert _python(Cm, Cm, 4) == ((False, None), 0)
    code, bound, cert = _c(Cm, Cm, 4)
    assert code == _tran_branch_c.QUIET and cert == 1
    assert bound == pytest.approx(1e-6 / np.sqrt(3), rel=1e-12)


@pytest.mark.parametrize('why, Cm, C0', [
    ## the rank holds (every singular value 1.5 tol) but 1/||B^-1||_F is
    ## 1.5 tol / sqrt(3), under tol
    ('bound below tol', _blocks(_sv([1.5e-9] * 3)), None),
    ## an entry between the two blocks the generic pattern made
    ('entry between blocks', None, None),
])
def test_what_the_certificate_cannot_prove_the_svd_decides(why, Cm, C0):
    if why == 'entry between blocks':
        C0 = _blocks(1e-6 * _Q)
        Cm = C0.copy()
        Cm[1, 3] = 1e-3
    C0 = Cm if C0 is None else C0
    verdict, svds = _python(Cm, C0, 4)
    assert verdict == (False, None) and svds == 1
    code, _bound, cert = _c(Cm, C0, 4)
    assert code == _tran_branch_c.SVD and cert == 0


def test_a_rank_that_drops_is_never_proved():
    """One singular value under tol: the SVD fires, with the null direction;
    neither certificate proves it (it cannot: its bound is under the
    smallest singular value)."""
    Cm = _blocks(_sv([1e-6, 1e-6, 0.5e-9]))
    verdict, svds = _python(Cm, Cm, 4)
    assert verdict[0] is True and verdict[1] is not None and svds == 1
    code, _bound, cert = _c(Cm, Cm, 4)
    assert code == _tran_branch_c.SVD and cert == 0
    tol = 1e-9
    data = _Check(Cm, 4)._branch_cert_data()
    assert not _certify(Cm, data, tol, _BranchCheck.BRANCH_CERT_GROWTH)[0]


def test_a_singular_submatrix_is_no_proof():
    """A 2 x 2 block with two equal rows (exactly singular: dgesv's info,
    numpy's LinAlgError), r0 counted as 4 by hand: the SVD decides, and
    fires, on both sides."""
    B3 = np.array([[1.0, 2.0, 0.0], [1.0, 2.0, 0.0], [0.0, 0.0, 1.0]]) * 1e-7
    Cm = _blocks(B3)
    verdict, svds = _python(Cm, Cm, 4)
    assert svds == 1 and verdict[0] is True
    assert _c(Cm, Cm, 4)[0] == _tran_branch_c.SVD


def _drawn(rng, sizes, u):
    """C with a 1 x 1 block 1.0 (the maximum) and blocks of `sizes`, each a
    random matrix scaled so its bound is tol (1 + u): the margins about the
    threshold the certificate decides at."""
    n = 2 + sum(sizes)
    Cm = np.zeros((n, n))
    Cm[1, 1] = 1.0
    o = 2
    for k in sizes:
        M = rng.standard_normal((k, k)) + np.eye(k) * 0.1
        b0 = 1.0 / np.linalg.norm(np.linalg.inv(M)) if k > 1 else abs(M[0, 0])
        Cm[o:o + k, o:o + k] = M * (1e-9 * (1.0 + u) / b0)
        o += k
    ## (every block's entries far under 1e-6 of the maximum: the cheap
    ## level fires)
    return Cm


@pytest.mark.parametrize('threads', [None, 1])
def test_the_c_certificate_is_the_python_one_bit_for_bit(threads):
    """Blocks of 1 to 150 unknowns (across OpenBLAS's threading thresholds),
    each at a margin about tol(1 + delta) from far under to far over: the C
    proves exactly what the Python proves, its smallest bound the same
    double; on the core's pattern (any size) and read whole (small)."""
    rng = np.random.default_rng(20261014)
    seen = {True: 0, False: 0}
    with threadpool_limits(threads):
        for sizes in ([1, 2], [3, 5, 8], [17], [33, 1, 2], [64], [99], [100], [101, 3], [128],
                      [150]):
            for u in (-0.5, -1e-6, -1e-12, 0.0, 1e-12, 1e-9, 1e-6, 1e-3, 0.5, 10.0):
                Cm = _drawn(rng, sizes, u)
                r0 = Cm.shape[0] - 1
                data = _Check(Cm, r0)._branch_cert_data()
                ok, worst = _certify(Cm, data, 1e-9, _BranchCheck.BRANCH_CERT_GROWTH)
                pat = np.flatnonzero(Cm).astype(np.int64)
                runs = [_c(Cm, Cm, r0, pat)] + ([_c(Cm, Cm, r0)] if Cm.size <= 4096 else [])
                for code, bound, cert in runs:
                    assert (code == _tran_branch_c.QUIET) == ok, (sizes, u)
                    assert cert == int(ok)
                    if ok:
                        assert float(bound).hex() == float(worst).hex(), (sizes, u)
                seen[ok] += 1
    assert seen[True] > 20 and seen[False] > 20, seen


## -- whole runs ---------------------------------------------------------------

def _rc():
    """Two nodes: 1 pF at `a`, 0.1 aF at `b` and between them -- the cheap
    level fires on every step (b's diagonal 2e-7 of a's), the rank holds:
    the certificate proves it, on its 2 x 2 block."""
    c = SubCircuit()
    c['vin'] = VPulse('in', gnd, v1=0.0, v2=1.0, td=1e-9, tr=1e-9, tf=1e-9, pw=5e-9, per=2e-8)
    c['r1'] = R('in', 'a', r=1e3)
    c['c1'] = C('a', gnd, c=1e-12)
    c['r2'] = R('a', 'b', r=1e3)
    c['c2'] = C('b', gnd, c=1e-19)
    c['c3'] = C('a', 'b', c=1e-19)
    return c


_PATHS = {'python': (False, False), 'python loop, C screen': (True, False), 'C loop': (True, True)}


def _run(monkeypatch, path, cert, build=_rc, **kw):
    screen, loop = _PATHS[path]
    monkeypatch.setattr(_tran_branch_c, 'ENABLED', screen)
    monkeypatch.setattr(_tran_loop_c, 'ENABLED', loop)
    before = _paths.snapshot()
    with warnings.catch_warnings():
        warnings.simplefilter('ignore')
        tr = Transient(build(), toolkit=circuit.numeric)
        tr.BRANCH_CERT = cert
        res = tr.solve(tend=3e-8, timestep=2e-10, **kw)
    st = tr.statistics
    return tr, (np.asarray(res.x, float).tobytes(), st.accepted_steps, st.rejected_steps,
                st.newton_iterations, getattr(st, 'branch_screens', 0)), _paths.since(before)


def test_a_run_is_the_same_with_the_certificate_on_every_path(monkeypatch):
    """The circuit the certificate proves: the waveform and the counts the
    same with it off and on, on the three paths; on, every screen the C
    makes is proved, and the C loop takes the attempts it handed back."""
    _tr, off, off_counts = _run(monkeypatch, 'C loop', False)
    assert off_counts.get('loop_c:bail:screen', 0) > 0
    for path in _PATHS:
        _tr, on, counts = _run(monkeypatch, path, True)
        assert on == off, path
        if path != 'python':
            assert counts['branch_c:cert'] == counts['branch_c:served'] > 0, (path, counts)
    _tr, _on, counts = _run(monkeypatch, 'C loop', True)
    assert counts.get('loop_c:bail:screen', 0) == 0
    assert counts['loop_c:attempts'] > 0.9 * off[1]


def test_its_data_are_made_once_and_kept_with_the_structural_rank(monkeypatch):
    """`_cert_data` at the first screen the cheap level fires on, kept --
    keyed by the structural rank's tuple, which the transient keeps between
    solves -- and made again where that rank is made again.  (This circuit's
    C is constant: after its first screen the screen's constant-matrix
    verdict answers, `_branch_const`.)"""
    calls = []
    real = _tran_branch._cert_data

    def counted(*a):
        calls.append(1)
        return real(*a)

    monkeypatch.setattr(_tran_branch, '_cert_data', counted)
    tr, _on, _counts = _run(monkeypatch, 'python', True)
    assert len(calls) == 1 and tr._branch_cert[1] is not None
    with warnings.catch_warnings():
        warnings.simplefilter('ignore')
        tr.solve(tend=1e-8, timestep=2e-10)
    assert len(calls) == 1
    del tr._branch_rank0
    tr._branch_structural_rank()
    assert tr._branch_cert_data() is not None and len(calls) == 2


def test_a_gold_deck_runs_in_the_c_loop(monkeypatch):
    """mosrect (the 27's; its cheap level fired on every step): the C loop
    takes its attempts with the certificate, the run the same as without."""
    from pycircuit.circuit import spice_import
    deck = benchdata.spice_data('Netlists/CircuitSim90/MOS2/mosrect.cir')
    if deck is None:
        pytest.skip('benchmark data not fetched (benchmarks/fetch_spice_suite.py)')

    def run(cert):
        circuit.default_toolkit = circuit.numeric
        imp = spice_import.import_netlist(deck)
        tr, skw = imp.transient()
        tr.BRANCH_CERT = cert
        before = _paths.snapshot()
        with warnings.catch_warnings():
            warnings.simplefilter('ignore')
            res = tr.solve(**skw)
        return np.asarray(res.x, float).tobytes(), _paths.since(before)

    off, off_counts = run(False)
    on, counts = run(True)
    assert on == off
    assert off_counts.get('loop_c:attempts', 0) == 0
    assert counts['loop_c:attempts'] > 900 and counts['branch_c:cert'] > 900

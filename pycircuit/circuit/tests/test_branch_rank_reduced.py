"""THE STRUCTURAL RANK OF C WITHOUT THE REFERENCE NODE (`_tran_branch.
_reduced_count`, `_BranchCheck._branch_structural_rank`): at every size the
branch check counts the singular values of C with the reference node's row
and column out -- the matrix the step solves -- above BRANCH_SCREEN_TOL of
its largest entry off the reference's diagonal, by blocks, with no band.
The count is numpy's dense SVD's of that reduced C at the run's own three
generic states, below and past `_packed.DENSE_MAX`; where it parts from the
full matrix's (a singular value between the two), the reduced one."""
import warnings

import numpy as np
import pytest

from pycircuit.circuit import _paths, circuit
from pycircuit.circuit.elements import IS, VS, C, L, R, SubCircuit, gnd
from pycircuit.circuit.tests.test_first_attempt_c import _chain
from pycircuit.circuit.tests.test_klu_order import _layers
from pycircuit.circuit.transient import Transient


def _coupled(k=12):
    """A ladder with capacitors between nodes as well as to ground, and an
    inductor (a branch row C holds): blocks larger than one."""
    c = SubCircuit()
    c['vs'] = VS('n0', gnd, v=1.0)
    for i in range(k):
        c[f'r{i}'] = R(f'n{i}', f'n{i + 1}', r=1e3)
        c[f'cg{i}'] = C(f'n{i + 1}', gnd, c=1e-12 * (1 + i % 3))
        if i % 2:
            c[f'cc{i}'] = C(f'n{i}', f'n{i + 1}', c=3e-13)
    c['l'] = L(f'n{k}', 'x', L=1e-9)
    c['rl'] = R('x', gnd, r=50.0)
    c['is'] = IS(f'n{k // 2}', gnd, i=1e-3)
    return c


def _numpy_rank0(tr):
    """`(r0, scale)` as numpy's dense SVD makes them: at the structural
    rank's three generic states, each C's largest entry off the reference's
    diagonal, the count of the reduced C's singular values above
    BRANCH_SCREEN_TOL of it, the largest of each."""
    rng = np.random.RandomState(20260910)
    g = tr.irefnode
    best, scale = 0, 0.0
    for _ in range(3):
        Cd = np.asarray(tr._C_fresh(rng.uniform(-1.0, 1.0, tr.cir.n)), dtype=float)
        off = np.abs(Cd)
        off[g, g] = 0.0
        nrm = float(off.max())
        keep = np.r_[0:g, g + 1:Cd.shape[0]]
        sv = np.linalg.svd(Cd[np.ix_(keep, keep)], compute_uv=False)
        scale, best = max(scale, nrm), max(best, int((sv > tr.BRANCH_SCREEN_TOL * nrm).sum()))
    return best, scale


@pytest.mark.parametrize('build', [lambda: _layers(8), lambda: _layers(16), lambda: _chain(400),
                                   _coupled], ids=['layers8', 'layers16', 'chain400', 'coupled'])
@pytest.mark.parametrize('cert', [None, 0], ids=['cert', 'held'])
def test_the_structural_rank_is_the_reduced_cs_count(build, cert):
    """The run's r0 and scale numpy's of the reduced C (`_numpy_rank0`),
    each state counted `branch.rank:reduced` -- with the certificate's
    generic C kept, and past its size (the states on one held array)."""
    with warnings.catch_warnings():
        warnings.simplefilter('ignore')
        tr = Transient(build(), toolkit=circuit.numeric)
    if cert is not None:
        tr.BRANCH_CERT_MAX_N = cert
    before = _paths.snapshot()
    got = tr._branch_structural_rank()
    counts = _paths.since(before)
    assert got == _numpy_rank0(tr) and got[0] > 0
    assert counts.get('branch.rank:reduced', 0) == 3 and not [k for k in counts if k in (
        'branch.rank:bounded', 'branch.rank:dense')], counts


def _run_rank(monkeypatch, tol_rel, dense_max):
    """`_coupled`'s run on THE PACKED LAYOUT at `BRANCH_SCREEN_TOL` `tol_rel`
    and `_packed.DENSE_MAX` `dense_max`: `(r0, the branch.rank counts, the
    check's own warnings, its error)`."""
    from pycircuit.circuit import _packed
    monkeypatch.setattr(_packed, 'PACKED_FROM', 0)
    monkeypatch.setattr(_packed, 'DENSE_MAX', dense_max)
    before = _paths.snapshot()
    with warnings.catch_warnings(record=True) as W:
        warnings.simplefilter('always')
        tr = Transient(_coupled(), toolkit=circuit.numeric)
        tr.BRANCH_SCREEN_TOL = tol_rel
        if dense_max <= tr.BRANCH_CERT_MAX_N:
            ## (as at size: the certificate's data are dense, kept below
            ## BRANCH_CERT_MAX_N, itself below DENSE_MAX)
            tr.BRANCH_CERT_MAX_N = 0
        tr.solve(tend=5e-10, timestep=1e-10, fixed_timestep=True)
    counts = {k: v for k, v in _paths.since(before).items() if k.startswith('branch.rank')}
    off = [str(w.message) for w in W if 'branch check' in str(w.message)]
    ## (no rank where the check switched itself off: the warning says why)
    return tr.__dict__.get('_branch_rank0', (None,))[0], counts, off, getattr(tr, '_branch_error', None)


@pytest.mark.parametrize('dense_max', [1 << 14, 1], ids=['below-dense-max', 'past-dense-max'])
def test_the_count_is_the_reduced_cs_at_every_size(dense_max, monkeypatch):
    """Below and past `_packed.DENSE_MAX` (here from one unknown) alike: at
    a tolerance a whole matrix's rounding margin `8 n eps ||C||` would pass
    (1e-15 of the largest entry at 17 unknowns -- as at 563 000 unknowns
    and the default) the check on, r0 the reduced matrix's dense count; and
    at a tolerance between the reduced and the full matrices' second
    singular values (the full counting 2, the reduced 1), the reduced 1."""
    with warnings.catch_warnings():
        warnings.simplefilter('ignore')
        tr = Transient(_coupled(), toolkit=circuit.numeric)
    Cd = np.asarray(tr._C_fresh(np.zeros(tr.cir.n)), dtype=float)
    g = tr.irefnode
    keep = np.r_[0:g, g + 1:Cd.shape[0]]
    off = Cd.copy()
    off[g, g] = 0.0
    scale = float(np.abs(off).max())
    full = np.linalg.svd(Cd, compute_uv=False)
    reduced = np.linalg.svd(Cd[np.ix_(keep, keep)], compute_uv=False)
    tiny = 1e-15
    assert 8 * Cd.shape[0] * np.finfo(float).eps * np.linalg.norm(Cd) > tiny * scale
    want = int((reduced > tiny * scale).sum())
    r0, counts, warned, err = _run_rank(monkeypatch, tiny, dense_max)
    assert (r0, counts.get('branch.rank:reduced'), warned, err) == (want, 3, [], None), (r0, counts, warned, err)
    ## (1e-2 of the largest entry lies between the reduced and the full
    ## matrices' second singular values)
    tol = 1e-2 * scale
    assert int((full > tol).sum()) == 2 and int((reduced > tol).sum()) == 1
    r0, counts, warned, err = _run_rank(monkeypatch, 1e-2, dense_max)
    assert (r0, counts.get('branch.rank:reduced'), warned, err) == (1, 3, [], None), (r0, counts, warned, err)


def test_a_non_finite_c_raises_and_the_check_says_so(monkeypatch):
    """A state's C not finite: `_reduced_count` raises (where a run's
    screen asked, the check then switches itself off with its warning) --
    not a silent count."""
    from pycircuit.circuit import _tran_branch
    with warnings.catch_warnings():
        warnings.simplefilter('ignore')
        tr = Transient(_coupled(), toolkit=circuit.numeric)
    real = tr._C_fresh

    def poisoned(x):
        Cd = np.array(real(x), dtype=float)
        Cd[0, 0] = np.nan
        return Cd
    monkeypatch.setattr(tr, '_C_fresh', poisoned)
    with pytest.raises(FloatingPointError, match='not finite'):
        tr._branch_structural_rank()
    nz = (np.array([0]), np.array([0]), np.array([np.inf]))
    with pytest.raises(FloatingPointError):
        _tran_branch._reduced_count(np.eye(1), nz, 1e-9, -1)

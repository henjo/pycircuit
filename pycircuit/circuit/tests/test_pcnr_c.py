"""PCNR's DC solve in C (`_pcnr_c`, `csrc/pcnrdc.c`) against its Python
reference, `pcnr.solve_dc`: the same solve with the C on and off, bit for
bit -- the answer, the limited unknowns and the iterations, or the
exception and its text -- and the C having SERVED it (a comparison the C
declined would compare the Python with itself).  The gold decks' census,
every rung of every deck, is the PCNR-in-C plan's C2 (`doc/transient_
runtime.md`); these are the circuits the suite carries."""
import itertools
import warnings

import numpy as np
import pytest
import sympy

from pycircuit._testing import klumem
from pycircuit.circuit import _paths, _pcnr_c, circuit, gnd, pcnr
from pycircuit.circuit import elements_hdl as eh
from pycircuit.circuit.dcanalysis import DC
from pycircuit.circuit.elements import VS, Diode, R, SubCircuit
from pycircuit.circuit.hdl import Behavioural, Branch, Contribution, limit_identity, var
from pycircuit.circuit.linearsolver import AutoSolver, DenseSolver
from pycircuit.circuit.tests.test_dc_pcnr import _bjt_mirror_20v
from pycircuit.utilities.param import Parameter

#: a level-3 card (`test_mos_gate_charge`'s) with its bulk junctions on
CARD = {'vto': 0.62, 'nsub': 1.066e16, 'kp': 6.3e-5, 'gamma': 0.64, 'phi': 0.31,
        'tox': 225e-10, 'vmax': 177269.0, 'l': 2e-6, 'w': 10e-6, 'IS': 1e-14}


def _cascode(vdd, vg2, vg1, drive_g2=True):
    """Two level-3 MOSFETs in a cascode, the bulk below ground: SPICE's
    four probes a device, the redundant `(b,d)` decided by `limit_block`'s
    Kruskal rule.  Without `drive_g2` the upper gate floats: the reduced
    matrix is singular."""
    c = SubCircuit()
    c['vdd'] = VS('vdd', gnd, v=vdd)
    if drive_g2:
        c['vg2'] = VS('g2', gnd, v=vg2)
    c['vg1'] = VS('g1', gnd, v=vg1)
    c['vb'] = VS('bulk', gnd, v=-0.5)
    c['rl'] = R('vdd', 'out', r=10e3)
    c['M2'] = eh.MosLevel3Hdl('out', 'g2', 'mid', 'bulk', **CARD)
    c['M1'] = eh.MosLevel3Hdl('mid', 'g1', gnd, 'bulk', **CARD)
    return c


def _ladder(stages=250, r1=100.0):
    """A source, then `stages` x (a series resistor, a SPICE diode to
    ground): 252 unknowns, so the DC's AutoSolver takes StructuralKLU.
    `r1` the first resistor."""
    c = SubCircuit()
    c['v'] = VS('n0', gnd, v=5.0)
    for k in range(1, stages + 1):
        c[f'r{k}'] = R(f'n{k - 1}', f'n{k}', r=r1 if k == 1 else 100.0)
        c[f'd{k}'] = eh.DiodeSpiceHdl(f'n{k}', gnd, IS=1e-14 * (1 + k % 7))
    return c


def _start(c, v):
    x0 = np.full(c.n, float(v))
    x0[c.get_node_index(gnd)] = 0.0
    return x0


def _solve(c, on, monkeypatch, **kw):
    """`pcnr.solve_dc` with the C on or off: `(the outcome, the C's
    counts during it)` -- the outcome `('ok', x bytes, v_lim bytes,
    iterations)` or `('raise', type, text)`."""
    monkeypatch.setattr(_pcnr_c, 'ENABLED', on)
    before = dict(_paths.COUNTS)
    kw.setdefault('linsolver', AutoSolver())
    with warnings.catch_warnings():
        warnings.simplefilter('ignore')
        try:
            x, v, its = pcnr.solve_dc(c, gnd, **kw)
            out = ('ok', x.tobytes(), v.tobytes(), its)
        except Exception as e:                           # noqa: BLE001
            out = ('raise', type(e), str(e))
    counts = {k: v - before.get(k, 0) for k, v in _paths.COUNTS.items()
              if k.startswith('pcnr_c:') and v != before.get(k, 0)}
    return out, counts


def _twin(c, monkeypatch, **kw):
    """The solve C on and C off: the C's outcome, asserted equal to the
    Python's and served."""
    got, counts = _solve(c, True, monkeypatch, **kw)
    ref, _ = _solve(c, False, monkeypatch, **kw)
    assert got == ref, (got[:1] + got[3:], ref[:1] + ref[3:])
    assert counts.get('pcnr_c:served') == 1, counts
    return got, counts


## ------------------------------------------------------------------------
## The refine kernel.

@pytest.mark.parametrize('which', ['mos3', 'gummel-poon'])
def test_the_refine_kernel_is_the_devices_limit(which):
    """Each class's refine kernel (`_pcnr_c.render_refine`) gives
    `PcnrDevice.limit` bit for bit: the law parameters (von at the last
    accepted sub-vector), the laws, and for a redundant probe the Kruskal
    and potentials of `limit_block` -- on draws where the laws bite."""
    c = _cascode(5.0, 3.0, 1.5) if which == 'mos3' else _bjt_mirror_20v()[0]
    dev = pcnr.pcnr_devices(c)[0]
    params, epar = dev.params(), circuit.defaultepar
    rng = np.random.default_rng(3)
    scale = np.array([0.8 if k == 'pnj' else 5.0 for k in dev.kinds])
    bit = 0
    for _ in range(500):
        vo = rng.uniform(-1, 1, dev.m) * scale
        vn = vo + rng.standard_normal(dev.m) * scale * 10.0 ** rng.uniform(-3, 1.5)
        xs = rng.uniform(-6, 6, dev.t)
        ref = dev.limit(vn.copy(), vo.copy(), params, epar, xs.copy())
        got = _pcnr_c.refine(dev, vn, vo, epar, xs)
        assert got is not None
        assert np.array_equal(np.asarray(ref).view(np.int64), got.view(np.int64))
        bit += not np.array_equal(ref, vn)
    assert bit > 100, bit


def test_the_refine_chains_from_the_cache_are_the_printed_ones(tmp_path, monkeypatch):
    """`render_refine`'s parameter chains -- the SymPy printing -- are kept
    beside the class's compiled info (`_pcnr_c._chains`): the first render
    prints them, the second reads them back without printing, and both
    are the text rendered with the cache off, for each class here."""
    from pycircuit.circuit import _hdl_cache, _hdl_climit
    monkeypatch.setenv('PYCIRCUIT_HDL_CACHE_DIR', str(tmp_path))
    monkeypatch.setattr(_hdl_cache, 'ENABLED', True)
    real = _hdl_climit._parameter_chains
    for c in (_cascode(5.0, 3.0, 1.5), _bjt_mirror_20v()[0]):
        cls = type(pcnr.pcnr_devices(c)[0].element)
        assert cls.__dict__.get('_hdl_cache_key'), cls
        printed = []

        def spy(*a, _printed=printed):
            _printed.append(1)
            return real(*a)
        monkeypatch.setattr(_hdl_climit, '_parameter_chains', spy)
        first, second = _pcnr_c.render_refine(cls), _pcnr_c.render_refine(cls)
        assert len(printed) == 1
        monkeypatch.setattr(_hdl_climit, '_parameter_chains', real)
        monkeypatch.setattr(_hdl_cache, 'ENABLED', False)
        ref = _pcnr_c.render_refine(cls)
        monkeypatch.setattr(_hdl_cache, 'ENABLED', True)
        assert first == second == ref


## ------------------------------------------------------------------------
## The solve.

@pytest.mark.parametrize('vdd,vg2,vg1,start', [
    (vdd, vg2, vg1, s) for (vdd, vg2, vg1), s in itertools.product(
        [(5.0, 3.0, 1.5), (20.0, 6.0, 0.5), (5.0, 1.0, 2.0)], [0.0, 20.0, -3.0])])
def test_a_cascode_solve_in_c_is_the_pythons(vdd, vg2, vg1, start, monkeypatch):
    _twin(_cascode(vdd, vg2, vg1), monkeypatch, x0=_start(_cascode(vdd, vg2, vg1), start))


@pytest.mark.parametrize('rev', [False, True])
@pytest.mark.parametrize('deform', ['plain', 'gshunt', 'source', 'maxiter'])
def test_a_mirror_solve_in_c_is_the_pythons(rev, deform, monkeypatch):
    """The BJT mirror from its 20 V start, both instance orders, under each
    deformation `solve_dc_continued` makes -- and the iteration limit,
    served as the Python's RuntimeError, text and all."""
    c, x0 = _bjt_mirror_20v(rev)
    kw = {'plain': {}, 'gshunt': {'gshunt': 1e-3}, 'source': {'source': 0.3, 'x0': None},
          'maxiter': {'maxiter': 3}}[deform]
    kw.setdefault('x0', x0)
    got, counts = _twin(c, monkeypatch, **kw)
    if deform == 'maxiter':
        assert got[0] == 'raise' and got[1] is RuntimeError and counts['pcnr_c:maxiter'] == 1
    else:
        assert got[0] == 'ok'


def test_a_singular_solve_is_numpys_error(monkeypatch):
    """A floating gate: the reduced matrix is singular at the first solve
    and numpy's LinAlgError is served, its own text."""
    c = _cascode(5.0, 3.0, 1.5, drive_g2=False)
    got, counts = _twin(c, monkeypatch, x0=_start(c, 0.0))
    assert got[0] == 'raise' and got[1] is np.linalg.LinAlgError, got
    assert counts['pcnr_c:singular'] == 1


def test_a_non_finite_iterate_is_the_pythons_error(monkeypatch):
    """A start with a NaN where no device probes (the supply node): the
    assembly and the dense solve carry it quietly, as numpy does, x_new is
    not finite, and the Python's RuntimeError (`pcnr.NONFINITE`) is served
    at the first iteration -- not handed back."""
    c = _cascode(5.0, 3.0, 1.5)
    x0 = _start(c, 1.0)
    x0[c.get_node_index('vdd')] = np.nan
    got, counts = _twin(c, monkeypatch, x0=x0)
    assert got == ('raise', RuntimeError, pcnr.NONFINITE % 1), got
    assert counts['pcnr_c:nonfinite'] == 1


class _Root(Behavioural):
    """A conductance through a square root: NaN at a negative probe -- a
    kernel's NaN at a finite point, as a wild iterate's (`exp` overflowing
    into `inf * 0`)."""

    @staticmethod
    def analog(p, n):
        b = Branch(p, n)
        r = var(sympy.sqrt(limit_identity(b.V)), 'r')
        return (Contribution(b.I, 1e-3 * r),)


def _ladder_and_root():
    c = _ladder()
    c['q'] = _Root('z', 'n100')
    return c


def _klu(on, monkeypatch, build, x0, first=None):
    """`_solve` on a fresh DC's own solver -- the AutoSolver over the
    pattern, the DC held, since the solver reaches its pattern through it
    weakly -- after a solve from `first` where one is given (KLU's numeric
    then held: the next solve refactors).  `(outcome, the C's counts,
    StructuralKLU's counters, whether it holds a numeric and the bytes KLU's
    objects hold: a numeric the C replaced and did not free shows there
    alone)`."""
    from pycircuit.circuit.dcanalysis import DC
    c = build()
    dc = DC(c, toolkit=circuit.numeric, pcnr=True)
    ls = dc._get_linearsolver()
    if first is not None:
        assert _solve(c, on, monkeypatch, x0=first, linsolver=ls)[0][0] == 'ok'
    got, counts = _solve(c, on, monkeypatch, x0=x0, linsolver=ls)
    sk = ls._choice
    assert type(sk).__name__ == 'StructuralKLU'
    return got, counts, (sk.refactors, sk.factors, sk.pivot_fallbacks, sk.residual_fallbacks,
                         sk._numeric is None, klumem.memusage(sk))


def _poisoned(build, node, value, rest=0.5):
    c = build()
    x0 = _start(c, rest)
    x0[c.get_node_index(node)] = value
    return x0


@pytest.mark.parametrize('case', ['rhs', 'matrix', 'refactor'])
def test_a_nan_in_the_klu_solve_is_the_pythons_error(case, monkeypatch):
    """On StructuralKLU a NaN goes through KLU's sequence in C -- the same
    library on the same values, its residual guard propagating the NaN as
    numpy's max does -- and the Python's RuntimeError (`pcnr.NONFINITE`) is
    SERVED: a NaN on the rhs (a NaN start on the source node), in the matrix
    (`_Root` at a negative probe), and on a held numeric (the refactor, then
    the residual's fresh factor), StructuralKLU's counters the same.  Handed
    back until 2026-10-19: the Meyer mode's e1480 reran five failing rungs a
    run in the Python's loop (0.3 .. 0.8 s each; three ended in KLU's failed
    factor, served the same way, `Singular matrix`)."""
    if case == 'matrix':
        build, x0, first = _ladder_and_root, _poisoned(_ladder_and_root, 'z', -0.5), None
    else:
        build, x0 = _ladder, _poisoned(_ladder, 'n0', np.nan)
        first = _start(_ladder(), 0.5) if case == 'refactor' else None
    got, counts, klu = _klu(True, monkeypatch, build, x0, first)
    ref, _counts, klu_ref = _klu(False, monkeypatch, build, x0, first)
    assert got == ref == ('raise', RuntimeError, pcnr.NONFINITE % 1), (got, ref)
    assert klu == klu_ref, (klu, klu_ref)
    assert counts.get('pcnr_c:nonfinite') == 1 and 'pcnr_c:bail:flags' not in counts, counts
    if case == 'refactor':
        assert klu[0] >= 1 and klu[3] >= 1, klu          # (refactored, then refactored again)


def test_an_infinite_klu_rhs_is_the_pythons_error(monkeypatch):
    """An infinity goes through KLU's sequence in C as a NaN does -- its
    residual's dense test hands back where the Python's `arr @ x - rhs`
    would warn -- and the Python's RuntimeError is SERVED, StructuralKLU's
    counters the same.  Handed back until 2026-11-08 (unmeasured, the
    residual's possible warning; roadmap 298)."""
    x0 = _poisoned(_ladder, 'n0', np.inf)
    got, counts, klu = _klu(True, monkeypatch, _ladder, x0)
    ref, _counts, klu_ref = _klu(False, monkeypatch, _ladder, x0)
    assert got == ref == ('raise', RuntimeError, pcnr.NONFINITE % 1), (got, ref)
    assert klu == klu_ref, (klu, klu_ref)
    assert counts.get('pcnr_c:nonfinite') == 1 and 'pcnr_c:bail:flags' not in counts, counts


def _active_chain(stages=800):
    """A pulsed inverter chain whose stages gain above one (vto 0.5, kp
    2e-4): its DC fails PCNR's plain solve, and on the ladders' rungs an
    iterate runs away -- a device derivative of ~1e92 times a limited
    residual of ~1e255, the Schur sum's rhs correction overflowing."""
    from pycircuit.circuit.elements import C, VPulse
    c = SubCircuit()
    c.add_node('vdd')
    c['vdd'] = VS('vdd', gnd, v=1.8)
    c.add_node('g0')
    c['vin'] = VPulse('g0', gnd, v1=0.0, v2=1.8, td=5e-9, tr=2e-9, tf=2e-9, pw=4e-8, per=1e-7)
    for k in range(stages):
        c.add_node(f'd{k}')
        c[f'rl{k}'] = R('vdd', f'd{k}', r=2e4)
        c[f'M{k}'] = eh.MosLevel1Hdl(f'd{k}', 'g0' if k == 0 else f'd{k - 1}', gnd, gnd,
                                     vto=0.5, kp=2e-4)
        c[f'c{k}'] = C(f'd{k}', gnd, c=2e-14)
    return c


def test_a_schur_overflow_at_a_run_away_rung_is_the_pythons_error(monkeypatch):
    """The active chain's DC, its PCNR calls recorded: the first one the C
    ends in KLU's failed factor, replayed on a fresh circuit with the C on
    and off -- the Python's sum overflows unwarned (`schur_reduce`'s
    errstate; warnings raised as errors here), the C carries the infinity on
    and serves its `Singular matrix`, StructuralKLU's counters the same.
    Handed back on the sum's overflow until 2026-11-08 (roadmap 298): the
    Python's dense loop reran it to the same failure, 2.3 s a rung at
    8 000 stages."""
    calls = []
    real = pcnr.solve_dc

    def recorded(*a, **k):
        before = dict(_paths.COUNTS)
        seed = (a, {n: (v.copy() if isinstance(v, np.ndarray) else v) for n, v in k.items()})
        try:
            return real(*a, **k)
        except Exception:
            if _paths.COUNTS.get('pcnr_c:singular', 0) > before.get('pcnr_c:singular', 0):
                calls.append(seed)
            raise
    monkeypatch.setattr(pcnr, 'solve_dc', recorded)
    with warnings.catch_warnings():
        warnings.simplefilter('ignore')
        DC(_active_chain(), toolkit=circuit.numeric).solve()
    monkeypatch.setattr(pcnr, 'solve_dc', real)
    assert calls, 'no PCNR call ended in a failed factor'
    (args, kw) = calls[0]
    sums = []
    schur = pcnr.schur_reduce

    def watched(*a, **k):
        f, J = schur(*a, **k)
        sums.append(bool(np.isinf(f).any()))
        return f, J
    out = {}
    for on in (True, False):
        c = _active_chain()
        kk = dict(kw, linsolver=DC(c, toolkit=circuit.numeric)._get_linearsolver())
        monkeypatch.setattr(_pcnr_c, 'ENABLED', on)
        if not on:
            monkeypatch.setattr(pcnr, 'schur_reduce', watched)
        before = dict(_paths.COUNTS)
        with warnings.catch_warnings():
            warnings.simplefilter('error', RuntimeWarning)
            try:
                pcnr.solve_dc(c, *args[1:], **kk)
                got = ('ok',)
            except Exception as e:                           # noqa: BLE001
                got = ('raise', type(e), str(e))
        sk = kk['linsolver']._choice
        out[on] = (got, (sk.refactors, sk.factors, sk.pivot_fallbacks, sk.residual_fallbacks),
                   {n: v - before.get(n, 0) for n, v in _paths.COUNTS.items()
                    if n.startswith('pcnr_c:') and v != before.get(n, 0)})
    assert out[True][:2] == out[False][:2], out
    assert out[True][0] == ('raise', np.linalg.LinAlgError, 'Singular matrix'), out
    assert any(sums), 'the Python\'s Schur sum did not overflow'
    counts = out[True][2]
    assert counts.get('pcnr_c:singular') == 1 and 'pcnr_c:bail:flags' not in counts, counts


@pytest.mark.parametrize('pattern', [True, False], ids=['pattern', 'no-pattern'])
def test_the_klu_solve_in_c_is_the_pythons(pattern, monkeypatch):
    """Above the AutoSolver's floor PCNR solves on StructuralKLU: the C
    stops at its first Schur matrix for the AutoSolver's own choice
    (`_choose`), then runs KLU's sequence -- on the core's pattern, or
    without one (`_tran_core.PATTERN` off: R made dense from J, the
    residual's dense test on it as built).  Two fresh DC analyses, C on
    and off: every PCNR call's outcome and StructuralKLU's counters and
    numeric after it the same."""
    from pycircuit.circuit import _tran_core
    from pycircuit.circuit.dcanalysis import DC
    monkeypatch.setattr(_tran_core, 'PATTERN', pattern)

    def run(on):
        _pcnr_c.ENABLED = on
        c = _ladder()
        dc = DC(c, toolkit=circuit.numeric, pcnr=True)
        log, real = [], pcnr.solve_dc

        def spy(*a, **k):
            try:
                x, v, its = real(*a, **k)
                log.append(('ok', x.tobytes(), v.tobytes(), its))
            except Exception as e:
                log.append(('raise', type(e), str(e)))
                raise
            sk = dc._default_linearsolver._choice
            log.append((type(sk).__name__, sk.refactors, sk.factors, sk.pivot_fallbacks,
                        sk.residual_fallbacks, sk._numeric is None, sk.analyses, klumem.memusage(sk)))
            return x, v, its
        pcnr.solve_dc = spy
        try:
            before = _paths.COUNTS['pcnr_c:served']
            with warnings.catch_warnings():
                warnings.simplefilter('ignore')
                x = np.asarray(dc.solve().x, dtype=float)
            rec = c.__dict__.get('_pcnr_c')
            return (log, x.tobytes(), _paths.COUNTS['pcnr_c:served'] - before,
                    rec is not None and rec.P is not None and bool(rec.kmaps))
        finally:
            pcnr.solve_dc = real
            _pcnr_c.ENABLED = True
    got, ref = run(True), run(False)
    assert got[0] == ref[0] and got[1] == ref[1]
    assert got[0][1][0] == 'StructuralKLU', got[0][1]
    assert got[2] == len(got[0]) // 2 and ref[2] == 0
    ## (the C on the pattern: KLU's entries gathered where they stand, the
    ## verdicts the sparse sums' -- pcnrdc.h's THE PATTERN)
    assert got[3] is pattern, 'the C did not solve on the pattern' if pattern else 'a pattern'


def test_a_structural_floor_below_the_sparse_one_is_served(monkeypatch):
    """An AutoSolver whose structural floor is below `min_n` (a DC asking
    for KLU from 40 unknowns; a 62-unknown ladder between the floors) has
    not chosen at the first solve, and the choice reads the matrix: the C
    stops at the first Schur matrix for it ('choose') and serves the solve,
    the same as the Python's.  Until 2026-10-19 the C made the choice on a
    zero matrix (fill 0: an unanalysed StructuralKLU), declined it
    ('linsolver'), and the Python then solved on a choice its own first
    matrix had not made."""
    from pycircuit.circuit.dcanalysis import DC
    monkeypatch.setattr(DC, '_STRUCTURAL_KLU_MIN_N', 40)
    out = []
    for on in (True, False):
        c = _ladder(60)
        dc = DC(c, toolkit=circuit.numeric, pcnr=True)
        ls = dc._get_linearsolver()
        assert ls.structural_min_n == 40 and c.n - 1 < ls.min_n
        got, counts = _solve(c, on, monkeypatch, x0=_start(c, 0.5), linsolver=ls)
        out.append((got, counts, type(ls._choice).__name__))
    (got, counts, ch), (ref, _c, ch_ref) = out
    assert got == ref and got[0] == 'ok', (got[:1], ref[:1])
    assert ch == ch_ref == 'StructuralKLU', (ch, ch_ref)
    assert counts.get('pcnr_c:served') == 1, counts


def test_the_first_reduced_matrix_holds_the_schur_entries_off_the_cores_pattern(monkeypatch):
    """At size the AutoSolver's choice reads the first reduced Schur matrix
    on pattern pages (`_reduced_J`): the dense one's bytes, `J[np.ix_(keep,
    keep)]`, the Schur entries off the core's pattern among them -- on a
    MOS chain a drain row's entry at its gate's column (the view's core
    holds no element across them), there not zero.  Until 2026-11-09 the
    matrix was made from the core's pattern alone, short of them."""
    from pycircuit.circuit import _tran_core
    from pycircuit.circuit.dcanalysis import DC
    monkeypatch.setattr(_tran_core, 'PATTERN_PAGES_FROM', 1)
    monkeypatch.setattr(DC, '_STRUCTURAL_KLU_MIN_N', 40)
    seen, real = [], _pcnr_c._reduced_J

    def spy(rec):
        R = real(rec)
        n, P, keep = rec.n, rec.P, rec.red >= 0
        off = ~np.isin(P, np.asarray(rec.core.arrays['pat'])) & keep[P // n] & keep[P % n]
        seen.append((R, _pcnr_c._dense_J(rec)[np.ix_(keep, keep)],
                     np.count_nonzero(rec.buf['J'][:len(P)][off])))
        return R
    monkeypatch.setattr(_pcnr_c, '_reduced_J', spy)
    c = _active_chain(60)
    ls = DC(c, toolkit=circuit.numeric, pcnr=True)._get_linearsolver()
    _out, counts = _solve(c, True, monkeypatch, x0=_start(c, 0.9), linsolver=ls)
    assert len(seen) == 1, counts
    R, D, nz_off = seen[0]
    assert nz_off > 0, 'no Schur entry off the core\'s pattern: the test reads none'
    assert _tran_core.pattern_of(R) is not None
    assert R.tobytes() == D.tobytes()
    assert type(ls._choice).__name__ == 'StructuralKLU'


def test_klus_entries_are_read_at_their_positions_in_packed_J(monkeypatch):
    """With a pattern J is packed on P (pcnrdc.h's THE PATTERN): each KLU
    entry is read at its position in it -- one off P at J's last, the zero
    -- and the foreign test at the positions of P's entries off KLU's
    pattern.  The circuits' KLU patterns hold P and no more, so a take made
    for the test: KLU's own short of a few entries (foreign, then) and with
    a few of the reduced matrix's entries off P."""
    from pycircuit.circuit.dcanalysis import DC
    monkeypatch.setattr(DC, '_STRUCTURAL_KLU_MIN_N', 40)
    c = _ladder(60)
    ls = DC(c, toolkit=circuit.numeric, pcnr=True)._get_linearsolver()
    _out, counts = _solve(c, True, monkeypatch, x0=_start(c, 0.5), linsolver=ls)
    rec = c.__dict__['_pcnr_c']
    assert counts.get('pcnr_c:served') == 1 and rec.P is not None and rec.klu is not None
    n, m, P, red = rec.n, int(rec.s.m), rec.P, rec.red
    unred = np.flatnonzero(red >= 0)
    take = rec.klu_keep[1]

    def full(t):
        return unred[t // m] * n + unred[t % m]
    on = np.isin(full(np.arange(m * m)), P)
    off_p = np.flatnonzero(~on)[:5]
    fake = np.sort(np.concatenate([take[3:], off_p]))
    rec.klu_keep = (rec.klu_keep[0], fake, *rec.klu_keep[2:])
    rec.kmaps.clear()
    rec.on_pattern(())
    _ch, kfull, ffull = rec.kmaps[(id(rec.klu), ())][:3]
    want = {f: i for i, f in enumerate(P.tolist())}
    assert kfull.tolist() == [want.get(f, len(P)) for f in full(fake).tolist()]
    assert np.count_nonzero(kfull == len(P)) == len(off_p) == 5
    kept = (red[P // n] >= 0) & (red[P % n] >= 0)
    foreign = P[kept][~np.isin(red[P[kept] // n] * m + red[P[kept] % n], fake)]
    assert len(foreign) >= 3 and P[ffull].tolist() == foreign.tolist()
    assert int(rec.s.nfor) == len(foreign)


## ------------------------------------------------------------------------
## On pattern pages: the C's n x n matrices resident where written.

def _stiff_ladder(stages=250):
    """`_ladder` with a first resistor of 1e-300 ohm: that row's sum
    |J| |x| past 1e300, where neither the convergence test's sparse verdict
    nor the residual's decides -- both dense tests made, solve after solve
    (the DC then fails: the point found does not solve its equations)."""
    return _ladder(stages, r1=1e-300)


def _dc_calls(build, on, monkeypatch):
    """A fresh DC on `build()`, the C on or off: each `pcnr.solve_dc` call's
    outcome and StructuralKLU's counters after it, the DC's own outcome, the
    C's counts and the circuit."""
    from pycircuit.circuit.dcanalysis import DC
    monkeypatch.setattr(_pcnr_c, 'ENABLED', on)
    c = build()
    dc = DC(c, toolkit=circuit.numeric, pcnr=True)
    log, real = [], pcnr.solve_dc

    def spy(*a, **k):
        try:
            x, v, its = real(*a, **k)
            log.append(('ok', x.tobytes(), v.tobytes(), its))
        except Exception as e:
            log.append(('raise', type(e), str(e)))
            raise
        finally:
            sk = dc._default_linearsolver._choice
            log.append((type(sk).__name__, getattr(sk, 'refactors', None),
                        getattr(sk, 'factors', None), getattr(sk, 'residual_fallbacks', None)))
        return x, v, its
    monkeypatch.setattr(pcnr, 'solve_dc', spy)
    before = dict(_paths.COUNTS)
    with warnings.catch_warnings():
        warnings.simplefilter('ignore')
        try:
            log.append(('dc', np.asarray(dc.solve().x, dtype=float).tobytes()))
        except Exception as e:                           # noqa: BLE001
            log.append(('dc raise', type(e)))
    monkeypatch.setattr(pcnr, 'solve_dc', real)
    counts = {k: v - before.get(k, 0) for k, v in _paths.COUNTS.items()
              if k.startswith('pcnr_c:') and v != before.get(k, 0)}
    return log, counts, c


def _pcnr_mats(c):
    """The circuit's PCNR record's n x n matrices: |J|, R and its core's G
    bin (J is packed on the pattern: pcnrdc.h's THE PATTERN)."""
    rec = c.__dict__['_pcnr_c']
    mats = {k: rec.buf[k] for k in ('Jabs', 'R')}
    mats['binG'] = rec.core.arrays['binG']
    return rec, mats


def test_on_pattern_pages_the_dense_tests_are_the_pythons(monkeypatch):
    """On pattern pages (`_tran_core.paged`, here from any size) the C's
    |J| and R and its core's G bin are pattern pages written on the pattern
    alone, J packed on it: R held zero -- the residual's dense test scatters
    KLU's entries in and takes them out after (pcx_kluseq.h's PCX_A_ZERO)
    -- and |J| written on the pattern where the convergence test goes
    dense.  On `_stiff_ladder`, both dense tests made solve after solve,
    every call is the Python's: its outcome and StructuralKLU's counters the
    same with the C on and off, each call served."""
    from pycircuit.circuit import _tran_core
    from pycircuit.circuit.dcanalysis import DC
    from pycircuit.circuit.tests.test_radau_pattern_pages import _paged
    monkeypatch.setattr(_tran_core, 'PATTERN_PAGES_FROM', 1)
    ## (StructuralKLU from 40 unknowns: the C's KLU on the pattern)
    monkeypatch.setattr(DC, '_STRUCTURAL_KLU_MIN_N', 40)
    got, counts, c = _dc_calls(lambda: _stiff_ladder(60), True, monkeypatch)
    ref, _counts, _c = _dc_calls(lambda: _stiff_ladder(60), False, monkeypatch)
    assert got[1][0] == 'StructuralKLU', got[1]
    assert got == ref
    calls = sum(1 for e in got if e[0] in ('ok', 'raise'))
    assert calls > 1 and counts.get('pcnr_c:served') == calls, (calls, counts)
    assert counts.get('pcnr_c:dense_conv', 0) > 0 and counts.get('pcnr_c:dense_resid', 0) > 0
    rec, mats = _pcnr_mats(c)
    assert all(_paged(v) for v in mats.values()), {k: _paged(v) for k, v in mats.items()}
    off = np.ones(c.n * c.n, dtype=bool)
    off[rec.P] = False
    assert np.count_nonzero(mats['Jabs']) > 0 and not np.any(mats['Jabs'][off])
    assert not np.any(mats['R'])
    J = rec.buf['J']
    assert J.shape == (len(rec.P) + 1,) and J[-1] == 0.0 and np.count_nonzero(J) > 0


def _resident_share(a):
    """The share of `a`'s pages resident and the process's own -- written:
    a read of a page never written maps the shared zero page
    (/proc/self/pagemap: present, bit 63; exclusively mapped, bit 56) --
    or None where the platform has no such map."""
    ptr = a.__array_interface__['data'][0]
    lo, hi = ptr // 4096, (ptr + a.nbytes + 4095) // 4096
    try:
        with open('/proc/self/pagemap', 'rb', buffering=0) as f:
            f.seek(lo * 8)
            raw = np.frombuffer(f.read((hi - lo) * 8), dtype=np.uint64)
    except OSError:
        return None
    want = (np.uint64(1) << np.uint64(63)) | (np.uint64(1) << np.uint64(56))
    return np.count_nonzero((raw & want) == want) / (hi - lo)


def test_at_size_the_dcs_matrices_are_resident_where_written(monkeypatch):
    """At size (`_tran_core.paged`: private 4 KB pages, resident where
    written) the C's |J|, R and its core's G bin hold a page or two a row
    (J packed: none): on `_stiff_ladder` of 1 500 stages (rows of three pages), both
    dense tests made, each under 3/4 of its pages resident (0.34-0.58
    measured) -- where a write of every entry (a dense test's zeroed R or
    |J| made whole) or one entry under a huge page makes all of it."""
    from pycircuit.circuit import _tran_core
    _log, counts, c = _dc_calls(lambda: _stiff_ladder(1500), True, monkeypatch)
    assert counts.get('pcnr_c:dense_conv', 0) > 0 and counts.get('pcnr_c:dense_resid', 0) > 0
    rec, mats = _pcnr_mats(c)
    assert _tran_core.paged(rec.core)
    share = {k: _resident_share(v) for k, v in mats.items()}
    if None in share.values():
        pytest.skip('no page map here')
    assert all(0.0 < s < 0.75 for s in share.values()), share


@pytest.mark.parametrize('threads', ['one', 'pool'])
def test_numpys_dgemv_rounds_a_row_once_a_nonzero_term(threads):
    """The premise of the C's sparse verdicts (pcx_kluseq.h's
    `pcx_klu_resid_sparse`, pcnrdc.c's `cf_sparse`): numpy's matmul -- its
    dgemv -- rounds a row at most once a NONZERO term (a zero product joins
    exactly, fused or not), so a row with k nonzeros is within 4 (k + 2) eps
    T of the same row summed in order, products rounded (T = sum |a| |x|).
    Drawn across the kernels' sizes (row counts off their widths, 2 .. 3001),
    signed rows cancelling to near zero and nonnegative ones, on one BLAS
    thread and on the pool."""
    import contextlib

    from pycircuit.circuit import _blas
    rng = np.random.default_rng(20261019)
    eps = np.finfo(float).eps
    hold = _blas._single_threaded_blas() if threads == 'one' else contextlib.nullcontext()
    worst, rows, differ = 0.0, 0, 0
    with hold:
        for n in (2, 3, 5, 17, 63, 64, 65, 250, 251, 999, 2333, 3001):
            for signed in (True, False):
                A = np.zeros((n, n))
                cols = []
                for i in range(n):
                    k = int(rng.integers(1, min(n, 12) + 1))
                    c = np.sort(rng.choice(n, size=k, replace=False))
                    v = rng.normal(size=k) * 10.0 ** rng.integers(-6, 7, size=k)
                    if signed and k > 1:
                        v[-1] = -float(np.sum(v[:-1])) * (1 + 1e-12 * rng.normal())
                    A[i, c] = v if signed else np.abs(v)
                    cols.append(c)
                x = rng.normal(size=n) * 10.0 ** rng.integers(-3, 4, size=n)
                if not signed:
                    x = np.abs(x)
                D = A @ x
                for i in range(n):
                    S = T = 0.0
                    for j in cols[i]:
                        S += A[i, j] * x[j]
                        T += abs(A[i, j]) * abs(x[j])
                    bound = 4.0 * (len(cols[i]) + 2) * eps * T
                    assert abs(D[i] - S) <= bound, (n, signed, i, D[i], S, T)
                    rows += 1
                    differ += D[i] != S
                    if T:
                        worst = max(worst, abs(D[i] - S) / bound)
    ## (alive: the dgemv sums in its own order -- measured 2026-10-19, 29 %
    ## of rows differ from the in-order sum, the worst 2.8 eps T, 0.07 of
    ## the bound)
    assert differ > rows // 10, (differ, rows)
    assert worst < 1.0, worst


## ------------------------------------------------------------------------
## The declines: counted, and the Python answers.

@pytest.mark.parametrize('why', ['gmin', 'underflow', 'scalar', 'linsolver', 'class', 'off'])
def test_each_decline_is_counted_and_the_python_answers(why, monkeypatch):
    from pycircuit.circuit import elements
    from pycircuit.circuit.tests.test_device_limiter import _cascode4, _mos4
    kw = {}
    c = _cascode(5.0, 3.0, 1.5)
    if why == 'gmin':
        kw['gmin'] = 1e-12
    elif why == 'scalar':
        c = SubCircuit()
        c['v'] = VS('a', gnd, v=1.0)
        c['r'] = R('a', 'b', r=1e3)
        c['d'] = elements.Diode('b', gnd)
    elif why == 'linsolver':
        kw['linsolver'] = DenseSolver()
    elif why == 'class':
        ## a model with no C kernels: pinned to numpy (a flat one has them
        ## since roadmap 168)
        from pycircuit.circuit import hdl
        cls = _mos4('group')
        hdl.set_backend('numpy', cls)
        c = _cascode4(cls, 5.0, 3.0, 1.5)
    monkeypatch.setattr(_pcnr_c, 'ENABLED', why != 'off')
    ref, _ = _solve(c, False, monkeypatch, **kw)
    monkeypatch.setattr(_pcnr_c, 'ENABLED', why != 'off')
    before = dict(_paths.COUNTS)
    with np.errstate(under='warn' if why == 'underflow' else 'ignore'), \
            warnings.catch_warnings():
        warnings.simplefilter('ignore')
        x, v, its = pcnr.solve_dc(c, gnd, linsolver=kw.pop('linsolver', AutoSolver()), **kw)
    assert ('ok', x.tobytes(), v.tobytes(), its) == ref
    assert _paths.COUNTS['pcnr_c:' + why] == before.get('pcnr_c:' + why, 0) + 1
    assert _paths.COUNTS['pcnr_c:served'] == before.get('pcnr_c:served', 0)


@pytest.mark.parametrize('late', ['own', 'toolkit', 'class'])
def test_a_later_devices_decline_is_the_solves(late, monkeypatch):
    """The record's checks are each device's in turn (`_record`: a class's
    kernels asked once a call): a decline on the 41st of 60 diodes of one
    class -- its instance's own `i`, a toolkit of another type -- or of a
    class only that diode has, is still the solve's: counted, not served,
    the Python's answer."""
    def build():
        c = _ladder(60)
        el = c['d41']
        if late == 'own':
            el.i = el.i
        elif late == 'toolkit':
            sub = type('_Tk', (type(el.toolkit),), {})
            tk = object.__new__(sub)
            tk.__dict__.update(el.toolkit.__dict__)
            el.toolkit = tk
        else:
            el.__class__ = type('_Late', (type(el),), {})
        return c
    real = _pcnr_c._class_ok
    monkeypatch.setattr(_pcnr_c, '_class_ok',
                        lambda cls, t: 'class' if cls.__name__ == '_Late' else real(cls, t))
    c, c_ref = build(), build()
    got, counts = _solve(c, True, monkeypatch, x0=_start(c, 0.5))
    ref, _c = _solve(c_ref, False, monkeypatch, x0=_start(c_ref, 0.5))
    assert got == ref and got[0] == 'ok', (got[:1], ref[:1])
    assert counts.get('pcnr_c:' + late) == 1 and 'pcnr_c:served' not in counts, counts


def test_a_one_ulp_fault_in_the_reference_is_seen(monkeypatch):
    """The comparison can see one ulp: the reference's device current moved
    by one ulp on every stamp, and the twin test's two answers differ."""
    c = _cascode(5.0, 3.0, 1.5)
    good, _ = _solve(c, True, monkeypatch, x0=_start(c, 20.0))
    real = pcnr.PcnrDevice._kernel_stamp

    def nudged(self, v, params, epar):
        cur, block = real(self, v, params, epar)
        cur[0] = np.nextafter(cur[0], np.inf)
        return cur, block
    monkeypatch.setattr(pcnr.PcnrDevice, '_kernel_stamp', nudged)
    bad, _ = _solve(c, False, monkeypatch, x0=_start(c, 20.0))
    assert good[0] == bad[0] == 'ok' and good != bad


## ------------------------------------------------------------------------
## The view after the circuit changes (2026-10-14).  `pcnr_view` was cached on
## the participants and the temperature alone, so a circuit edited between
## two solves kept the first one's view: an added element was missing from
## the second solve (its answer the OLD circuit's, status 'used'), an added
## unknown left the C's core on a view one size short of its buffers (ASan:
## a heap-buffer-overflow in core.c's groups_into, from pcnrdc.c), and a
## participant's parameter moved the device but not its remainder's stub.

def _clamp():
    """A source, a resistor, a SPICE diode to ground."""
    c = SubCircuit()
    c['v'] = VS('a', gnd, v=5.0)
    c['r'] = R('a', 'b', r=1e3)
    c['d'] = eh.DiodeSpiceHdl('b', gnd, IS=1e-14)
    return c


def _with_resistor(c):
    """A resistor across the diode: the same unknowns."""
    c['r2'] = R('b', gnd, r=1e3)
    return c


def _with_ammeter(c):
    """A resistor across the diode through a 0 V source: a node and a branch
    more."""
    c['vm'] = VS('b', 'c', v=0.0)
    c['r2'] = R('c', gnd, r=1e3)
    return c


@pytest.mark.parametrize('on', [True, False])
@pytest.mark.parametrize('add', [_with_resistor, _with_ammeter])
def test_an_element_added_between_two_solves_is_in_the_second(add, on, monkeypatch):
    """The second solve is a fresh circuit's, served by the C where it is on."""
    c = _clamp()
    first, _ = _solve(c, on, monkeypatch)
    assert first[0] == 'ok'
    add(c)
    got, counts = _solve(c, on, monkeypatch)
    ref, _ = _solve(add(_clamp()), on, monkeypatch)
    assert got[0] == 'ok' and got[1] != first[1][:len(got[1])]
    assert got == ref
    assert counts.get('pcnr_c:served') == (1 if on else None), counts


def _common_emitter(re):
    """A Gummel-Poon stage: `re` is in the device's affine remainder, which
    the view stamps from its own stub."""
    c = SubCircuit()
    c['vcc'] = VS('vcc', gnd, v=5.0)
    c['rc'] = R('vcc', 'c', r=1e3)
    c['vb'] = VS('vb', gnd, v=1.0)
    c['rb'] = R('vb', 'b', r=10e3)
    c['q'] = eh.GummelPoonNpn4Hdl('c', 'b', gnd, gnd, re=re)
    return c


@pytest.mark.parametrize('on', [True, False])
def test_a_participants_parameter_changed_between_two_solves_is_in_the_second(on, monkeypatch):
    """`re` 1 -> 200 ohm between two solves read v(c) -0.682 V -- the device at
    the new card, its stub at the old -- where a fresh circuit reads 4.230."""
    c = _common_emitter(1.0)
    first, _ = _solve(c, on, monkeypatch)
    assert first[0] == 'ok'
    c['q'].ipar.re = 200.0
    got, counts = _solve(c, on, monkeypatch)
    ref, _ = _solve(_common_emitter(200.0), on, monkeypatch)
    assert got[0] == 'ok' and got[1] != first[1]
    assert got == ref
    assert counts.get('pcnr_c:served') == (1 if on else None), counts


def test_an_unchanged_circuit_keeps_its_view():
    """Solves of an unchanged circuit share one view (the core built on it is
    the C's record)."""
    c = _clamp()
    devs = pcnr.pcnr_devices(c)
    with warnings.catch_warnings():
        warnings.simplefilter('ignore')
        pcnr.solve_dc(c, gnd)
        view = pcnr.pcnr_view(c, devs)
        pcnr.solve_dc(c, gnd)
    assert pcnr.pcnr_view(c, devs) is view


def _with_a_participant(vdd_expr, nested):
    """A parameterised circuit -- `vdd` its parameter -- whose source takes
    its value through an expression, beside a diode (a PCNR participant,
    so a DC builds the PCNR view).  `nested`: the source inside a
    subcircuit instance whose own parameter is that expression."""
    tk = circuit.numeric

    class Top(SubCircuit):
        instparams = [Parameter('vdd', default=0.0)]  # noqa: RUF012

    class Inner(SubCircuit):
        instparams = [Parameter('v', default=0.0)]  # noqa: RUF012
        terminals = ('p', 'n')

    top = Top(vdd=2.0, toolkit=tk)
    n1, n2 = top.add_nodes('1', '2')
    if nested:
        inner = Inner(n1, gnd, v=vdd_expr, toolkit=tk)
        inner['V'] = VS('p', 'n', v='v', toolkit=tk)
        top['X'] = inner
    else:
        top['V'] = VS(n1, gnd, v=vdd_expr, toolkit=tk)
    top['R'] = R(n1, n2, r=1e3, toolkit=tk)
    top['D'] = Diode(n2, gnd, toolkit=tk)
    return top, n1


@pytest.mark.parametrize('nested', [False, True], ids=['top', 'nested'])
def test_a_parameter_from_a_parent_s_expression_survives_the_pcnr_view(nested):
    """A DC on a circuit with a PCNR participant builds the PCNR view, which
    holds the circuit's own elements; a parameter given as an expression of
    a parent's keeps its value through it, and the answer is the circuit's
    (until 2026-10-22 the view's build re-resolved every element against
    the view's own, empty, parameters: the source's `v` reset to its
    default 0 V -- the answer 0 V, and the circuit left so)."""
    c, n1 = _with_a_participant('vdd', nested)
    ref, _n1 = _with_a_participant(2.0, nested)
    got = DC(c).solve()
    want = DC(ref).solve()
    assert float(got.v(n1, gnd)) == 2.0
    assert np.asarray(got.x).tobytes() == np.asarray(want.x).tobytes()
    v = c['X']['V'] if nested else c['V']
    assert v.iparv.v == 2.0


def test_the_device_list_is_kept_until_a_participant_or_the_topology_changes():
    """`pcnr_devices` keeps its records on the circuit under the view's key:
    a source's value written keeps them, a participant's parameter or a new
    element brings new ones -- a new list each call either way."""
    c = _clamp()
    first = pcnr.pcnr_devices(c)
    again = pcnr.pcnr_devices(c)
    assert again is not first and all(a is b for a, b in zip(first, again))
    c['v'].ipar.v = 4.0
    assert all(a is b for a, b in zip(first, pcnr.pcnr_devices(c)))
    c['d'].ipar.IS = 2e-14
    moved = pcnr.pcnr_devices(c)
    assert not any(a is b for a, b in zip(first, moved))
    _with_resistor(c)
    assert not any(a is b for a, b in zip(moved, pcnr.pcnr_devices(c)))


def test_the_seed_s_refine_through_the_c_kernel_is_the_python_one(monkeypatch):
    """`pcnr.refine` runs each device's C refine kernel where it has one --
    the C solve's own -- and gives the Python law's answer to the bit; with
    the C off, the Python law."""
    c = _clamp()
    devs = pcnr.pcnr_devices(c)
    m = sum(d.m for d in devs)
    rng = np.random.default_rng(5)
    served = []
    real = _pcnr_c.refine

    def spy(*a, **k):
        out = real(*a, **k)
        served.append(out is not None)
        return out
    monkeypatch.setattr(_pcnr_c, 'refine', spy)
    for _ in range(20):
        x = rng.uniform(-1.0, 6.0, c.n)
        v_old, v_new = rng.uniform(-1.0, 3.0, m), rng.uniform(-2.0, 8.0, m)
        monkeypatch.setattr(_pcnr_c, 'ENABLED', True)
        on = pcnr.refine(devs, v_old, v_new, x_old=x)
        monkeypatch.setattr(_pcnr_c, 'ENABLED', False)
        off = pcnr.refine(devs, v_old, v_new, x_old=x)
        assert on.tobytes() == off.tobytes()
    assert served and all(served)


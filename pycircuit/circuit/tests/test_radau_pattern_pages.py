"""PATTERN PAGES FOR RADAU'S RECORDS (plan 261's L6a): at size
(`_tran_core.paged`: a core with a pattern, `PATTERN_PAGES_FROM` entries)
the n x n matrices a Radau run keeps -- the C transform loop's hand-back
records (`_tran_tcloop_c._record`) and its `_Geq`, the passes' C and G
(`_tran_core.passes`), the transform record's end sets and the copies it
hands back -- go on pattern pages (`_tran_core.pattern_zeros`: private 4 KB
pages, resident where written), where under huge pages a matrix written on
its pattern alone was resident whole; and the operating point's PCNR
scratch kept on the circuit is handed back after a transient's op
(`_pcnr_c.release`).  Each is the same bytes: the run is the dense run, bit
for bit.  So too the coupled Newton in C at size (plan 261's L14c,
`_tran_radau_c._put_paged`): the fallback analyses the stage KLU and stays in
C, writes every assembly's C and G to one slot, and the memo takes the last
assembly's alone, on pattern pages."""
import mmap
import warnings

import numpy as np
import pytest

from pycircuit.circuit import (
    _paths,
    _pcnr_c,
    _tran_core,
    _tran_radau_c,
    _tran_tcloop_c,
    circuit,
)
from pycircuit.circuit import elements_hdl as eh
from pycircuit.circuit.elements import VS, C, R, SubCircuit, VPulse, gnd
from pycircuit.circuit.integrator import Gear2Integrator, RadauIIA3Integrator
from pycircuit.circuit.transient import Transient


def _chain(n):
    """A pulsed MOS chain with gate-drain capacitors, one three stages back
    at every third stage, and a junction's depletion charge at every fifth
    node (C moving with the state)."""
    c = SubCircuit()
    c.add_node('vdd')
    c['vdd'] = VS('vdd', gnd, v=1.8)
    c.add_node('g0')
    c['vin'] = VPulse('g0', gnd, v1=0.0, v2=1.8, td=2e-9, tr=1e-9, tf=1e-9, pw=2e-8, per=5e-8)
    for k in range(n):
        g = 'g0' if k == 0 else f'd{k - 1}'
        c.add_node(f'd{k}')
        c[f'rl{k}'] = R('vdd', f'd{k}', r=2e4)
        c[f'M{k}'] = eh.MosLevel1Hdl(f'd{k}', g, gnd, gnd)
        c[f'c{k}'] = C(f'd{k}', gnd, c=2e-14)
        c[f'cm{k}'] = C(f'd{k}', g, c=2e-15)
        if k >= 3 and k % 3 == 0:
            c[f'cl{k}'] = C(f'd{k}', f'd{k - 3}', c=1e-15)
        if k % 5 == 0:
            c[f'D{k}'] = eh.DiodeSpiceHdl(gnd, f'd{k}', cjo=2e-14, vj=0.7, m=0.5)
    return c


def _b(a):
    return np.asarray(a, dtype=float).tobytes()


def _memo(tr):
    """The device memo's records: each state's arrays' bytes, by generation."""
    memo = tr.__dict__.get('_dev_memo') or ({}, {})
    return tuple(sorted((k, tuple(sorted((m, _b(v)) for m, v in r.items()))) for k, r in g.items())
                 for g in memo)


def _paged(a):
    """Whether `a`'s data is on pattern pages (a private mapping of its own)."""
    b = a
    while isinstance(b, np.ndarray) and b.base is not None:
        b = b.base
    if isinstance(b, memoryview):
        b = b.obj
    return isinstance(b, mmap.mmap)


def _run(monkeypatch, pages, integ=RadauIIA3Integrator, n=300, cir=None, **kw):
    """The chain's run with the pattern pages from any size, or from none:
    its bytes and statistics, the state it leaves (its matrices' bytes and
    whether they are on pages), its paths (the release's own left out)."""
    monkeypatch.setattr(_tran_core, 'PATTERN_PAGES_FROM', 1 if pages else 1 << 62)
    before = _paths.snapshot()
    with warnings.catch_warnings():
        warnings.simplefilter('ignore')
        tr = Transient(cir if cir is not None else _chain(n), toolkit=circuit.numeric,
                       integrator=integ())
        res = tr.solve(**(kw or {'tend': 4e-8, 'timestep': 1e-9}))
    st = tr.statistics
    run = (_b(res.x), _b(res.sweep_values), st.accepted_steps, st.rejected_steps,
           st.newton_iterations)
    mats = {k: getattr(tr, k, None) for k in ('_Cmat', '_Geq')}
    left = (tuple((k, _b(v), np.shape(v)) for k, v in mats.items()), _memo(tr),
            _b(tr._q_cache[1]), _b(tr._iq))
    on = {k: _paged(v) for k, v in mats.items() if isinstance(v, np.ndarray)}
    paths = {k: v for k, v in _paths.since(before).items()
             if not k.startswith('once:') and k != 'pcnr_c:released'}
    return run, left, on, paths, tr


def _moved(a, b):
    """The path counts that move with the pattern pages, taken out: each
    fresh C the screen reads (`_C_fresh`) the core's pass where it was
    `cir.C`'s assembly (`plan.C`, its two batches), each screen served on
    its pattern where it declined `large`, each frozen transform's factors
    made at the structure (`radau.frozen:packed`) where some were folded --
    the rest left to compare."""
    a, b = dict(a), dict(b)
    fresh = a.get('plan.C:served', 0) - b.get('plan.C:served', 0)
    assert fresh >= 0
    assert a.get('batch.C:served', 0) - b.get('batch.C:served', 0) == 2 * fresh
    assert b.get('core.p:served', 0) - a.get('core.p:served', 0) == fresh
    big = a.pop('branch_c:large', 0) - b.pop('branch_c:large', 0)
    assert b.get('branch_c:pattern', 0) - a.get('branch_c:pattern', 0) == big
    assert b.get('branch_c:served', 0) - a.get('branch_c:served', 0) == big
    ## (the frozen factors made at the structure's positions: no fold)
    assert 'radau.frozen:packed' not in a and 'radau.frozen:fold' not in b
    assert b.get('radau.frozen:packed', 0) >= a.get('radau.frozen:fold', 0)
    for k in ('plan.C:served', 'batch.C:served', 'core.p:served', 'branch_c:served',
              'branch_c:pattern', 'radau.frozen:packed', 'radau.frozen:fold',
              'radau.frozen:fold.same'):
        a.pop(k, None)
        b.pop(k, None)
    return {k: (a.get(k), b.get(k)) for k in set(a) | set(b) if a.get(k) != b.get(k)}


@pytest.mark.parametrize('kw', [{}, {'tend': 2e-8, 'timestep': 2e-10, 'fixed_timestep': True}],
                         ids=['adaptive', 'fixed'])
def test_on_pattern_pages_the_c_loop_is_the_python_loop(monkeypatch, kw):
    """With the pattern pages, the C transform loop's run is the Python
    loop's -- bytes, state, and every path count the C reproduces for the
    attempts it makes (the screens on the pattern, the fresh C's the core's
    pass)."""
    from pycircuit.circuit import _tran_tcloop_c
    with monkeypatch.context() as m:
        m.setattr(_tran_tcloop_c, 'ENABLED', False)
        a = _run(m, True, **kw)
    b = _run(monkeypatch, True, **kw)
    assert b[3].get('tcloop_c:steps', 0) > 10 and 'tcloop_c:steps' not in a[3]
    loop = ('tcloop_c:', 'loop_c:', 'radau.fuse:handed', 'radau.fuse:unhanded', 'sources_c:')
    pa = {k: v for k, v in a[3].items() if not k.startswith(loop)}
    pb = {k: v for k, v in b[3].items() if not k.startswith(loop)}
    assert a[0] == b[0] and a[1] == b[1]
    assert {k: (pa.get(k), pb.get(k)) for k in set(pa) | set(pb) if pa.get(k) != pb.get(k)} == {}


@pytest.fixture(autouse=True)
def _warm(monkeypatch):
    """A first run for the process's builds made once."""
    with monkeypatch.context() as m:
        _run(m, False, n=20, tend=4e-9, timestep=1e-9)


@pytest.mark.parametrize('integ, kw', [
    (RadauIIA3Integrator, {}),
    (RadauIIA3Integrator, {'tend': 2e-8, 'timestep': 2e-10, 'fixed_timestep': True}),
    (Gear2Integrator, {})], ids=['radau', 'radau-fixed', 'gear2'])
def test_a_run_on_pattern_pages_is_the_dense_run(monkeypatch, integ, kw):
    """A 305-unknown chain: with the pattern pages its run, its statistics,
    every path count, the state it leaves -- `_Cmat`, `_Geq`, the device
    memo's every record, the cached q and iq -- are the dense run's, bytes
    and shapes; Radau's matrices left on pattern pages, the C loop serving."""
    a = _run(monkeypatch, False, integ, **kw)
    b = _run(monkeypatch, True, integ, **kw)
    assert a[0] == b[0] and a[1] == b[1]
    assert _moved(a[3], b[3]) == {}
    if integ is RadauIIA3Integrator:
        assert b[3].get('tcloop_c:steps', 0) > 10, b[3]
        assert b[2] == {'_Cmat': True, '_Geq': True} and not any(a[2].values())
        recs = [v for g in b[4].__dict__['_dev_memo'] for r in g.values() for k, v in r.items()
                if k in 'CG']
        assert recs and all(_paged(v) for v in recs)


def _fallback_run(monkeypatch, pages, maxiter, reltol, tend=4e-8):
    """The chain with the Newton's iterations held so the transform stalls,
    the C transform loop off: each stall the coupled Newton.  `(its bytes and
    statistics, its paths, the transient)`."""
    monkeypatch.setattr(_tran_tcloop_c, 'ENABLED', False)
    monkeypatch.setattr(_tran_core, 'PATTERN_PAGES_FROM', 1 if pages else 1 << 62)
    before = _paths.snapshot()
    with warnings.catch_warnings():
        warnings.simplefilter('ignore')
        tr = Transient(_chain(300), toolkit=circuit.numeric, integrator=RadauIIA3Integrator(),
                       reltol=reltol)
        tr.par.maxiter = maxiter
        res = tr.solve(tend=tend, timestep=1e-9)
    st = tr.statistics
    return ((_b(res.x), _b(res.sweep_values), st.accepted_steps, st.rejected_steps),
            _paths.since(before), tr)


@pytest.mark.parametrize('maxiter, reltol', [(2, 3e-2), (3, 1e-3)])
def test_a_fallback_at_size_is_the_dense_run(monkeypatch, maxiter, reltol):
    """The coupled Newton at size, one slot and the last assembly's records:
    the run's bytes and statistics the dense run's, in C from the first
    fallback where the dense run's first is the Python's (`klu_unanalysed`)."""
    a, ca, _ta = _fallback_run(monkeypatch, False, maxiter, reltol)
    b, cb, _tb = _fallback_run(monkeypatch, True, maxiter, reltol)
    assert a == b
    assert cb.get('radau_c:served', 0) >= 1, cb
    assert ca.get('radau_c:klu_unanalysed') == 1 and not cb.get('radau_c:klu_unanalysed'), (ca, cb)
    assert cb['radau_c:served'] >= ca.get('radau_c:served', 0)


def test_a_fallback_at_size_keeps_the_last_assembly_s_matrices(monkeypatch):
    """Of a served call's records at size the last assembly's alone hold `C`
    and `G`, on pattern pages, the values a pass at its stages gives; every
    assembly's hold `q` and `i`; the slot and the core's G bin handed back,
    the C told one slot for that call alone."""
    calls = []
    real = _tran_radau_c._put_paged

    def spy(rec, put, o, A):
        assert A >= 2
        puts = []
        real(rec, lambda y, r: (puts.append((y.copy(), r)), put(y, r)), o, A)
        calls.append((A, puts))
        assert rec.s.one_slot == 0
        ## (the slot handed back: the middle stage's whole pages read zero)
        assert not np.count_nonzero(o['Ca'][0, 1]) and not np.count_nonzero(o['Ga'][0, 1])
    monkeypatch.setattr(_tran_radau_c, '_put_paged', spy)
    _, paths, tr = _fallback_run(monkeypatch, True, 3, 1e-3)
    assert calls, paths
    for A, puts in calls:
        assert len(puts) == 3 * A
        for y, r in puts[:3 * (A - 1)]:
            assert set(r) == {'q', 'i'}
        for y, r in puts[3 * (A - 1):]:
            assert set(r) == {'q', 'i', 'C', 'G'} and _paged(r['C']) and _paged(r['G'])
            assert _tran_core.pattern_of(r['C']) is not None
            P = _tran_core.passes(tr, y, 'CG')
            assert _b(r['C']) == _b(P['C']) and _b(r['G']) == _b(P['G'])


def test_the_passes_on_pattern_pages_are_the_dense_passes(monkeypatch):
    """`passes`: C and G on pattern pages at size, their bytes the dense
    call's; q and i the same."""
    _r, _l, _o, _p, tr = _run(monkeypatch, False, n=60, tend=2e-9, timestep=1e-9)
    x = np.asarray(tr._q_cache[0], dtype=float).copy()
    monkeypatch.setattr(_tran_core, 'PATTERN_PAGES_FROM', 1 << 62)
    a = _tran_core.passes(tr, x, 'GCiq')
    monkeypatch.setattr(_tran_core, 'PATTERN_PAGES_FROM', 1)
    b = _tran_core.passes(tr, x, 'GCiq')
    assert a is not None and b is not None
    assert all(_b(a[k]) == _b(b[k]) and a[k].shape == b[k].shape for k in 'GCiq')
    assert _paged(b['C']) and _paged(b['G']) and not _paged(a['C']) and not _paged(a['G'])
    assert np.count_nonzero(b['G']) > 0 and np.count_nonzero(b['C']) > 0


def test_the_ops_pcnr_scratch_is_handed_back_and_the_next_op_is_the_same(monkeypatch):
    """After a transient's operating point the circuit's PCNR record keeps
    its tables, its n x n scratch handed back (read as zeros: |J|, R, the
    core's C and G bins; J itself is packed on the pattern, written whole
    each assembly); a second transient on the same circuit -- its op on the
    handed-back record -- is the first one's."""
    monkeypatch.setattr(_tran_core, 'PATTERN_PAGES_FROM', 1)
    cir = _chain(80)
    before = _paths.snapshot()
    a = _run(monkeypatch, True, cir=cir, tend=1e-8, timestep=1e-9)
    rec = cir.__dict__.get('_pcnr_c')
    if not isinstance(rec, _pcnr_c._Rec):
        pytest.skip('the op is not the PCNR C record\'s')
    assert _paths.since(before).get('pcnr_c:released') == 1
    page = 4096
    assert rec.buf['J'].shape == (len(rec.P) + 1,)
    for k in ('Jabs', 'R', 'Cs', 'binG'):
        v = rec.core.arrays[k] if k == 'binG' else rec.buf[k]
        ptr = v.__array_interface__['data'][0]
        lo = (-ptr) % page // 8
        hi = (ptr + v.nbytes) // page * page - ptr
        if hi // 8 > lo:
            assert not np.any(v[lo:hi // 8]), k
    b = _run(monkeypatch, True, cir=cir, tend=1e-8, timestep=1e-9)
    assert a[0] == b[0] and a[1] == b[1]


def test_give_back_zeros_whole_pages_and_keeps_the_ends():
    """`give_back`: the whole pages of the data read as zeros, the partial
    pages at its ends kept."""
    page = 4096
    if getattr(mmap, 'MAP_PRIVATE', None) is None:
        pytest.skip('no private mappings here')
    ## (on a mapping of its own: nothing else on its end pages)
    m = mmap.mmap(-1, 8 * page, flags=mmap.MAP_PRIVATE | mmap.MAP_ANONYMOUS)
    a = np.frombuffer(m, dtype=np.float64, count=5 * page // 8 + 7, offset=24)
    a[:] = 1.5
    _tran_core.give_back(a)
    ptr = a.__array_interface__['data'][0]
    lo = (-ptr) % page // 8
    hi = ((ptr + a.nbytes) // page * page - ptr) // 8
    if _tran_core._libc() is None:
        pytest.skip('no madvise here')
    assert hi > lo
    assert not np.any(a[lo:hi]) and np.all(a[:lo] == 1.5) and np.all(a[hi:] == 1.5)


@pytest.mark.parametrize('c', [2.5, 0.0, -0.0, -1.5, np.inf, np.nan])
def test_a_scaled_matrix_on_pattern_pages_is_numpys(c):
    """`scaled_on_pattern`: `c * G`'s bytes -- on pattern pages where numpy's
    product is zero off the pattern as they are (`c` finite, its sign +),
    else numpy's own."""
    n = 64
    rng = np.random.default_rng(3)
    pat = np.flatnonzero(rng.random(n * n) < 0.05)
    G = np.zeros(n * n)
    G[pat] = rng.standard_normal(len(pat))
    G[pat[:3]] = -0.0
    G = G.reshape(n, n)
    with np.errstate(all='ignore'):
        want = c * G
        got = _tran_core.scaled_on_pattern(c, G, pat, n)
    assert _b(got) == _b(want) and got.shape == want.shape


def test_zero_pages_read_maps_no_memory():
    """The transform's `zero_pages`: private pages -- a read of all of them
    makes none resident (a shared mapping's read allocates each)."""
    from pycircuit.circuit import _tran_radau_tc
    if getattr(mmap, 'MAP_PRIVATE', None) is None:
        pytest.skip('no private mappings here')

    def rss():
        with open('/proc/self/statm') as f:
            return int(f.read().split()[1]) * 4096

    a = _tran_radau_tc.zero_pages(8 << 20)        # 64 MB
    r0 = rss()
    assert float(a.sum()) == 0.0
    assert rss() - r0 < 16 << 20


def _paged_matrix(n, rng, monkeypatch):
    """An n x n matrix on pattern pages of a random pattern -- entries in
    the branch rows and columns too, a NaN, signed zeros among them."""
    monkeypatch.setattr(_tran_core, 'PATTERN_PAGES_FROM', 1)
    pat = np.flatnonzero(rng.random(n * n) < 0.08).astype(np.int64)
    rng.shuffle(pat)
    full = np.zeros(n * n)
    full[pat] = rng.standard_normal(len(pat))
    full[pat[:3]] = 0.0
    full[pat[3]] = -0.0
    full[pat[4]] = np.nan
    return _tran_core.on_pattern(full.reshape(n, n), pat, n), full.reshape(n, n), pat


def test_counts_comparisons_and_couplings_on_pattern_pages_are_numpys(monkeypatch):
    """At a matrix's pattern: its nonzero count, its equality with another on
    the same pattern (and with one that is not), the node-to-node nonzeros
    `capacitive_loop_currents` joins, row by row -- numpy's on the whole."""
    from pycircuit.circuit import _tran_history
    rng = np.random.default_rng(11)
    n, nn = 60, 45
    A, full, pat = _paged_matrix(n, rng, monkeypatch)
    assert _tran_core.pattern_of(A) is pat
    assert _tran_core.pattern_count(A) == np.count_nonzero(full)
    B = _tran_core.on_pattern(np.nan_to_num(full), pat, n)
    D = _tran_core.on_pattern(np.where(np.isnan(full), 0.0, full), pat, n)
    for X, Y in ((B, D), (B, A), (A, A.copy())):
        assert _tran_core.pattern_equal(X, Y) == bool(np.array_equal(np.asarray(X), np.asarray(Y)))
    got = _tran_history._node_couplings(A, nn)
    want = np.nonzero(full[:nn, :nn] != 0.0)
    assert all(np.array_equal(g, w) for g, w in zip(got, want))


@pytest.mark.parametrize('extra', [False, True], ids=['copy', 'off-pattern'])
def test_the_loop_reads_a_memo_record_zero_off_the_pattern(monkeypatch, extra):
    """The C transform loop's entry takes x_n's C from the memo at the
    core's pattern only where the record is zero off it: a dense copy of
    each record's C serves as the pattern-page one does; with one nonzero
    off the pattern the entry declines (`tcloop_c:memo_pattern`)."""
    from pycircuit.circuit import _tran_tcloop_c
    monkeypatch.setattr(_tran_core, 'PATTERN_PAGES_FROM', 1)
    ## (the loop on its packed sets, as at size: the check is theirs)
    monkeypatch.setattr(_tran_tcloop_c, 'PACKED', 'always')
    r0 = _tran_tcloop_c.run
    done = []

    def run(loop):
        if not done:
            done.append(1)
            memo = loop.tr.__dict__.get('_dev_memo') or ({}, {})
            pat = set(np.asarray(_tran_core.core_for(loop.tr).arrays['pat']).tolist())
            for gen in memo:
                for rec in gen.values():
                    if 'C' in rec:
                        C = np.array(rec['C'], dtype=float)
                        if extra:
                            n = C.shape[0]
                            off = next(k for k in range(n * n) if k not in pat)
                            C.reshape(-1)[off] = 1e-30
                        rec['C'] = C
        return r0(loop)

    monkeypatch.setattr(_tran_tcloop_c, 'run', run)
    before = _paths.snapshot()
    with warnings.catch_warnings():
        warnings.simplefilter('ignore')
        tr = Transient(_chain(300), toolkit=circuit.numeric, integrator=RadauIIA3Integrator())
        tr.solve(tend=1e-8, timestep=1e-9)
    counts = _paths.since(before)
    assert done
    assert (counts.get('tcloop_c:memo_pattern', 0) > 0) == extra, counts
    if not extra:
        assert counts.get('tcloop_c:packed', 0) > 0 and counts.get('tcloop_c:steps', 0) > 0, counts


def test_a_matrix_on_pattern_pages_is_reduced_on_them(monkeypatch):
    """`remove_row_col` of a matrix on pattern pages -- the step's J the
    first attempts' error test reduces (plan 261's L19: the dense reduce
    wrote all (n-1)^2 entries of an IBM power grid's J, 5 GB resident): on
    pattern pages itself, its pattern the reduced one, the dense reduce's
    bytes -- the reference row first, inside and last; a NaN and signed
    zeros among the entries."""
    from pycircuit.circuit import _stamp_plan, analysis
    monkeypatch.setattr(_stamp_plan, 'PAGES_FROM', 1)
    rng = np.random.default_rng(5)
    n = 40
    A, full, _pat = _paged_matrix(n, rng, monkeypatch)
    for k in (0, 17, n - 1):
        got = analysis.remove_row_col((A,), k, circuit.numeric)[0]
        want = analysis.remove_row_col((full.copy(),), k, circuit.numeric)[0]
        assert _b(got) == _b(want), k
        assert _tran_core.pattern_of(got) is not None, k
        assert _tran_core.pattern_count(got) == np.count_nonzero(want), k


def test_the_plans_matrices_at_size_are_on_small_pages(monkeypatch):
    """The stamp plan's G and C at size (`_stamp_plan.PAGES_FROM` entries),
    and the DC's reduced G: each bin summed in bincount's order from +0.0
    onto small pages -- bincount's bytes, only the bins' pages written
    (bincount's own n^2 output, written where its bins land, was resident by
    huge pages: an IBM grid's C 2.7 GB) -- and no pattern kept with them (a
    caller may write them anywhere)."""
    from pycircuit.circuit import _stamp_plan
    cir = _chain(30)
    x = np.random.default_rng(3).uniform(-1.0, 1.0, cir.n)
    ep = circuit.defaultepar
    want = {m: _b(_stamp_plan.assemble_matrix(cir, m, x, (ep,))) for m in ('G', 'C')}
    wred = _b(_stamp_plan.assemble_matrix(cir, 'G', x, (ep,), drop=0))
    monkeypatch.setattr(_stamp_plan, 'PAGES_FROM', 1)
    monkeypatch.setattr(_tran_core, 'PATTERN_PAGES_FROM', 1)
    for m in ('G', 'C'):
        got = _stamp_plan.assemble_matrix(cir, m, x, (ep,))
        assert _b(got) == want[m], m
        assert _tran_core._mapping(got) is not None and _tran_core.pattern_of(got) is None, m
    red = _stamp_plan.assemble_matrix(cir, 'G', x, (ep,), drop=0)
    assert _b(red) == wred and _tran_core._mapping(red) is not None


def test_the_plans_pages_threshold_is_the_cores():
    """`_stamp_plan.PAGES_FROM`, read there without importing `_tran_core`
    (an assembly imports nothing), is `_tran_core.PATTERN_PAGES_FROM`."""
    from pycircuit.circuit import _stamp_plan
    assert _stamp_plan.PAGES_FROM == _tran_core.PATTERN_PAGES_FROM


def test_the_couplings_are_read_at_the_plans_positions(monkeypatch):
    """`capacitive_loop_currents` given the circuit, on a C its stamp plan
    made and not on pattern pages: the node-to-node nonzeros read at the
    plan's positions (`_tran_branch._plan_nonzeros`, proven all of them by
    their count) -- numpy's `nonzero` of the whole block, row by row, which
    at an IBM grid's 25 850 unknowns wrote a 0.6 GB mask."""
    from pycircuit.circuit import _tran_branch, _tran_history
    cir = _chain(30)
    x = np.random.default_rng(4).uniform(-1.0, 1.0, cir.n)
    C = np.asarray(cir.C(x))
    nn = cir.n - 3
    calls, real = [], _tran_branch._plan_nonzeros

    def counted(c, M):
        calls.append(1)
        return real(c, M)
    monkeypatch.setattr(_tran_branch, '_plan_nonzeros', counted)
    got = _tran_history._node_couplings(C, nn, cir)
    want = np.nonzero(C[:nn, :nn] != 0.0)
    assert calls and all(np.array_equal(g, w) for g, w in zip(got, want))


def test_the_newton_and_error_test_buffers_are_on_small_pages(monkeypatch):
    """At size the C Newton's record keeps its zeroed m x m buffers -- |J|,
    the sparse full Newton's R -- and the error test's R on small pages
    (`_tran_core.pattern_zeros`), each written on the pattern alone (lte.h's
    R HELD ZERO: its dense test zeroed all m^2 entries first): under huge
    pages each was resident whole, an IBM grid's run 6.5 GB through its
    loop."""
    monkeypatch.setattr(_tran_core, 'PATTERN_PAGES_FROM', 1)
    with warnings.catch_warnings():
        warnings.simplefilter('ignore')
        tr = Transient(_chain(60), toolkit=circuit.numeric, chord_jacobian=False)
        tr.solve(tend=2e-9, timestep=1e-10)
    rec = tr.__dict__['_newton_c']
    m = rec.n - 1
    assert rec.klu is not None
    assert _tran_core._mapping(rec.buf['absJ']) is not None
    big = [a for a in rec.klu_keep if type(a) is np.ndarray and a.size == m * m]
    assert big and all(_tran_core._mapping(a) is not None for a in big)
    ctx = tr.step_controller.__dict__.get('_lte_c')
    assert ctx is not None and ctx.klu_keep, 'the controller made no sparse error test'
    assert _tran_core._mapping(ctx.klu_keep[3]) is not None

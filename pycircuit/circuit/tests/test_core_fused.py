"""The evaluate core's groups: a class's kernel (`_hdl_cbackend.fuse_csrc`,
one masked kernel a class) called once a core call for every pass the call
wants, the class's batches scattering from its staging -- the answers,
statistics and warnings of the same core with every batch asking the
kernel for its own pass alone (the grouping patched off here).  A kernel
that is not the one the passes' `CKernel`s call, or none, is not grouped.
(The kernel's pass sets against each pass alone, on drawn and special
states: `test_twins_random.py`.)"""
import warnings

import numpy as np
import pytest

from pycircuit.circuit import _hdl_cbackend as cb
from pycircuit.circuit import circuit, hdl
from pycircuit.circuit import elements_hdl as eh
from pycircuit.circuit.elements import VS, R, gnd
from pycircuit.circuit.tests.test_hdl_cbackend import c_backend, needs_cc
from pycircuit.circuit.tests.test_newton_c import _gp_chain, _mos_chain
from pycircuit.circuit.tests.test_psp_limit_c import _chain, _stage
from pycircuit.circuit.transient import Transient

HEAD = 'void hdl_fn(const double *x, const double *p, double *out) {'


def _src(*body):
    return '\n'.join([HEAD, *body, '}']) + '\n'


def test_the_union_keeps_each_pass_and_guards_each_run():
    text, stats = cb.fuse_csrc({
        'G': _src('  const double L_a = x[0]*p[0];', '  const double L_b = L_a + 1.0;',
                  '  out[0] = L_b;', '  out[1] = L_a;'),
        'q': _src('  const double L_a = x[0]*p[0];', '  const double L_c = 2.0*L_a;',
                  '  out[0] = L_c;')})
    assert stats == {'union': 3, 'runs': 3, 'printed': 4}
    body = text.split('\n')
    assert body[0] == ('void hdl_fn(const double *x, const double *p, long want, '
                       'double *oG, double *oC, double *oi, double *oq) {')
    ## L_a in both passes (bits 1|8), L_b in G's alone, L_c in q's alone
    assert '  if (want & 9) {' in body and '    L_a = x[0]*p[0];' in body
    assert body.index('  if (want & 1) {') < body.index('    L_b = L_a + 1.0;')
    assert '    oG[1] = L_a;' in body and '    oq[0] = L_c;' in body


@pytest.mark.parametrize('srcs, why', [
    ({'G': _src('  const double L_a = x[0];', '  out[0] = L_a;'),
      'C': _src('  const double L_a = x[1];', '  out[0] = L_a;')}, 'two texts'),
    ({'G': _src('  const double L_a = x[0];', '  out[0] = L_a;'),
      'C': _src('  double L_a = x[0];', '  out[0] = L_a;')}, 'does not make'),
    ({'G': _src('  const double L_a = x[0];', '  out[1] = L_a;'),
      'C': _src('  const double L_a = x[0];', '  out[0] = L_a;')}, 'does not make'),
    ({'G': _src('  const double L_a = x[0];', '  out[0] = L_a;',
                '  const double L_b = x[1];'),
      'C': _src('  const double L_a = x[0];', '  out[0] = L_a;')}, 'does not make'),
    ({}, 'no pass'),
    ({'G': 'void f(double *x) {\n  out[0] = 1.0;\n}\n',
      'C': _src('  out[0] = 1.0;')}, 'printed frame'),
    ## a name read in a pass that does not print it
    ({'G': _src('  const double L_a = x[0];', '  const double L_b = 2.0*L_a;', '  out[0] = L_b;'),
      'C': _src('  const double L_b = 2.0*L_a;', '  out[0] = L_b;')}, 'outside its passes'),
    ({'G': _src('  const double L_a = x[0];', '  out[0] = L_a;'),
      'C': _src('  out[0] = L_a;')}, 'outside it'),
])
def test_the_union_refuses_what_it_cannot_print(srcs, why):
    with pytest.raises(cb.FuseRefused, match=why):
        cb.fuse_csrc(srcs)


def _no_groups(batches, b_pass):
    """`_tran_core._fuse_groups` grouping nothing: every batch asks its
    class's kernel for its own pass."""
    return [], [-1] * len(batches), [None] * len(batches)


def _run(build, fuse, monkeypatch, **kw):
    """One transient: grouped where `fuse`, else every batch a pass at a
    time (for this run only)."""
    from pycircuit.circuit import _tran_core
    make = kw.pop('make', {})
    with monkeypatch.context() as mp, warnings.catch_warnings(record=True) as W:
        if not fuse:
            mp.setattr(_tran_core, '_fuse_groups', _no_groups)
        warnings.simplefilter('always')
        ## (a finalizer's warning is not the run's)
        warnings.simplefilter('ignore', ResourceWarning)
        tr = Transient(build(), toolkit=circuit.numeric, **make)
        res = tr.solve(**kw)
    rec = tr.__dict__.get('_tran_core')
    core = rec[1] if rec is not None else None
    st_ = {k: getattr(tr.statistics, k) for k in tr.statistics.__slots__ if 'seconds' not in k}
    return (np.asarray(res.x, float).tobytes(), np.asarray(res.sweep_values, float).tobytes(),
            st_, sorted(str(w.message) for w in W), core, tr)


def _same(build, monkeypatch, **kw):
    a = _run(build, False, monkeypatch, **dict(kw))
    b = _run(build, True, monkeypatch, **dict(kw))
    assert a[0] == b[0] and a[1] == b[1], 'the solution moved'
    assert a[2] == b[2], ('the statistics moved', a[2], b[2])
    assert a[3] == b[3], 'the warnings moved'
    return a, b


def _mixed():
    """A MOS chain driving two Gummel-Poon stages: two fused groups."""
    c = _mos_chain(3)
    c.add_node('vcc')
    c['vcc'] = VS('vcc', gnd, v=3.0)
    prev = 'd2'
    for k in range(2):
        c.add_node(f'b{k}')
        c.add_node(f'c{k}')
        c[f'rb{k}'] = R(prev, f'b{k}', r=1e4)
        c[f'rc{k}'] = R('vcc', f'c{k}', r=1e3)
        c[f'Q{k}'] = eh.GummelPoonNpnHdl(f'c{k}', f'b{k}', gnd)
        prev = f'c{k}'
    return c


CASES = {
    'mos': (_mos_chain, {}),
    'gp-chord': (_gp_chain, {'chord_jacobian': True}),
    'psp-stage': (_stage, {}),
    'psp-chain-limited': (lambda: _chain(4, va=0.6), {}),
    'mixed': (_mixed, {}),
}


@needs_cc
@pytest.mark.parametrize('fixed', [True, False], ids=['fixed', 'adaptive'])
@pytest.mark.parametrize('case', list(CASES))
def test_a_grouped_core_gives_the_answers_of_a_pass_at_a_time(case, fixed, monkeypatch):
    build, make = CASES[case]
    kw = {'tend': 1e-6, 'timestep': 2e-8, 'make': make}
    if fixed:
        kw['fixed_timestep'] = True
    a, b = _same(build, monkeypatch, **kw)
    if b[4] is None:
        pytest.skip('the core does not serve this circuit here')
    assert b[4].nz >= 1 and all(b[4].z_fn[:b[4].nz]), 'a group not fused'
    assert a[4].nz == 0 or not any(a[4].z_fn[:a[4].nz])


def _alias(fk):
    """A `CKernel`'s function pointer to the same machine code as `fk`'s
    but another object: a kernel the core did not group with."""
    return fk.ffi.cast(fk.ffi.typeof(fk.cfn), fk.cfn)


@needs_cc
def test_a_stale_or_missing_fused_kernel_is_not_called(monkeypatch):
    """A pass's `CKernel` rebound behind the group's back (a pointer to the
    kernel's code, another object): the group's kernel is not the one
    `el.G` calls, and the batches call their own a pass at a time; the
    class's kernel gone from its info: the same."""
    a = _run(_mos_chain, False, monkeypatch, tend=4e-7, timestep=2e-8, fixed_timestep=True)
    b = _run(_mos_chain, True, monkeypatch, tend=4e-7, timestep=2e-8, fixed_timestep=True)
    core, tr = b[4], b[5]
    assert core is not None and core.nz == 1 and core.z_fn[0]
    info = core.zgroups[0][0]
    fn = info['funcs']['G']
    old = fn.__dict__['_hdl_c']
    other = cb.CKernel(old.ffi, _alias(info['_c_fused']), old.shape, old.nx,
                       (old.n_p, old.t_index), old.key, old.built_s, bit=old.bit)
    monkeypatch.setattr(fn, '_hdl_c', other)
    assert core.probe(tr.epar).__class__ is not str
    assert not core.z_fn[0], 'a group whose kernel el.G does not call was called'
    monkeypatch.setattr(fn, '_hdl_c', old)
    assert core.probe(tr.epar).__class__ is not str and core.z_fn[0]
    monkeypatch.setitem(info, '_c_fused', None)
    assert core.probe(tr.epar).__class__ is not str and not core.z_fn[0]
    monkeypatch.undo()
    ## and a transient on a core whose group fell back: the same bytes
    c = _run(_mos_chain, True, monkeypatch, tend=4e-7, timestep=2e-8, fixed_timestep=True)
    assert a[:4] == b[:4] == c[:4]


@needs_cc
def test_the_backend_switched_and_back_binds_a_new_fused_kernel(monkeypatch):
    """`set_backend` numpy then C: the class detaches (the core declines:
    unbound), then binds a new kernel and its passes' `CKernel`s, which the
    core takes after checking they call it."""
    cls = type(eh.MosLevel1Hdl('d', 'g', gnd, gnd))
    a = _run(_mos_chain, True, monkeypatch, tend=4e-7, timestep=2e-8, fixed_timestep=True)
    fk0 = cls._hdl_info['_c_fused']
    with c_backend(cls):
        hdl.set_backend('numpy', cls)
        assert cls._hdl_info['_c_fused'] is None
        hdl.set_backend('c', cls)
        fk1 = cls._hdl_info['_c_fused']
        assert fk1 is not None and fk1 is not fk0
        b = _run(_mos_chain, True, monkeypatch, tend=4e-7, timestep=2e-8, fixed_timestep=True)
        assert b[4].z_fn[0] and b[4].zk[0] is fk1
    assert a[:4] == b[:4]


@needs_cc
def test_a_class_without_its_fused_kernel_is_not_grouped(monkeypatch):
    """The grouping's refusal, counted: no kernel in the class's info (its
    passes' `CKernel`s still bound: the batches call them a pass at a
    time)."""
    from pycircuit.circuit import _paths
    cls = type(eh.MosLevel1Hdl('d', 'g', gnd, gnd))
    monkeypatch.setitem(cls._hdl_info, '_c_fused', None)
    monkeypatch.setitem(cls._hdl_info, '_c_fused_status', 'refused (a test)')
    before = _paths.snapshot()
    b = _run(_mos_chain, True, monkeypatch, tend=2e-7, timestep=2e-8, fixed_timestep=True)
    assert b[4].nz == 0
    assert _paths.since(before).get('once:core.fuse:nokernel', 0) >= 1


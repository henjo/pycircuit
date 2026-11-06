"""The evaluation core's packed sets (`PCX_WANT_SETS`, core.h's THE PACKED
SETS; plan 261's L4d2, roadmap 270): with the pattern, every matrix a call
reads or writes packed -- G's bins in `binGp`, C, the held Cin, Geq and J
the caller's, npat values in the pattern's order -- where a caller holding
no n x n set asks for it.  Each entry's operations are the pattern's call's
(the pull through `wpall`, each bin's index in the pattern; the companion
over the packed index): the packed values are the dense call's at the
pattern, bit for bit, and the dense bins are left as they were."""
import warnings

import numpy as np
import pytest

from pycircuit.circuit import _tran_core, circuit, gnd
from pycircuit.circuit import elements as E
from pycircuit.circuit import elements_hdl as eh
from pycircuit.circuit.circuit import SubCircuit
from pycircuit.circuit.transient import Transient

ABI = _tran_core._ABI


def _circuit(n=40):
    """A MOS chain with ground, gate-drain and long-range capacitors and a
    junction's depletion charge at every fifth node (C moving with x)."""
    c = SubCircuit()
    c['vdd'] = E.VS('vdd', gnd, v=1.8)
    c['vg'] = E.VPulse('g', gnd, v1=0.0, v2=1.8, td=2e-9, tr=1e-9, tf=1e-9, pw=5e-9, per=2e-8)
    for k in range(n):
        g = 'g' if k == 0 else f'd{k - 1}'
        c[f'R{k}'] = E.R('vdd', f'd{k}', r=5e3)
        c[f'C{k}'] = E.C(f'd{k}', gnd, c=2e-14)
        c[f'Cm{k}'] = E.C(f'd{k}', g, c=2e-15)
        c[f'M{k}'] = eh.MosLevel1Hdl(f'd{k}', g, gnd, gnd)
        if k >= 3 and k % 3 == 0:
            c[f'Cl{k}'] = E.C(f'd{k}', f'd{k - 3}', c=1e-15)
        if k % 5 == 0:
            c[f'D{k}'] = eh.DiodeSpiceHdl(gnd, f'd{k}', cjo=2e-14, vj=0.7, m=0.5)
    return c


@pytest.fixture(scope='module')
def core():
    with warnings.catch_warnings():
        warnings.simplefilter('ignore')
        tr = Transient(_circuit(), toolkit=circuit.numeric)
        tr.solve(tend=2e-9, timestep=1e-9)
    core = _tran_core.core_for(tr)
    if core is None or not core.npat:
        pytest.skip('the core does not serve this circuit on a pattern')
    return core, core.ready(tr.epar)


def _call(core, T, x, want, formula, sets, Cin, rng):
    """One call -- dense on the pattern, or packed -- its outputs and the
    dense G bin after it."""
    ffi, cfn, _ = _tran_core.driver()
    fb, dptr, n, npat = ffi.from_buffer, core.dptr, core.n, core.npat
    vec = [rng.standard_normal(n) for _ in range(4)]          # q1, q2, iq1, u
    size = npat if sets else n * n
    C, Geq, J = np.zeros(size), np.zeros(size), np.zeros(size)
    ## (the vectors from the same values in both calls: an entry a call does
    ## not write compares as itself)
    q, iq, F = np.full(n, 0.25), np.full(n, 0.5), np.full(n, 0.75)
    bits = want | ABI['PCX_WANT_PATTERN'] | (ABI['PCX_WANT_SETS'] if sets else 0)
    a = (0.7, -0.3, 0.05, 1.3e-10, 0.6)
    cfn(core.cs, fb(dptr, x), T, bits, formula, *a, *(fb(dptr, v) for v in vec),
        fb(dptr, Cin), fb(dptr, C), fb(dptr, q), fb(dptr, iq), fb(dptr, Geq), fb(dptr, F),
        fb(dptr, J))
    return C, Geq, J, q, iq, F


FORMULAS = ('PCX_FORMULA_GEAR2', 'PCX_FORMULA_EULER', 'PCX_FORMULA_TRAPEZOIDAL',
            'PCX_FORMULA_THETA', 'PCX_FORMULA_PASSES')
WANTS = {'all': 'PCX_WANT_ALL', 'G and I': ('PCX_PASS_G', 'PCX_PASS_I'),
         'C': ('PCX_PASS_C',), 'I': ('PCX_PASS_I',), 'G packed J': ('PCX_PASS_G', 'PCX_WANT_PACKED')}


def _want(key):
    names = WANTS[key]
    return ABI[names] if isinstance(names, str) else sum(ABI[k] for k in names)


@pytest.mark.parametrize('formula', FORMULAS)
@pytest.mark.parametrize('wkey', sorted(WANTS))
def test_the_packed_sets_are_the_pattern_calls_values(core, formula, wkey):
    """At random states, for every formula and want: C, Geq and J packed are
    the dense pattern call's at the pattern -- the held Cin packed where C
    is not made, J packed with or without Geq -- q, iq and F the same; G is
    `binGp`, the dense call's G bin at the pattern; the dense bin untouched
    by a packed call."""
    core, T = core
    pat, n = core.arrays['pat'], core.n
    rng = np.random.default_rng(hash((formula, wkey)) % 2**32)
    want = _want(wkey)
    for _ in range(3):
        x = rng.uniform(-1.0, 1.0, n)
        x[0] = 0.0
        Cin = np.zeros(n * n)
        Cin[pat] = rng.standard_normal(len(pat)) * 1e-14
        st = np.random.default_rng(7)
        d = _call(core, T, x, want, ABI[formula], False, Cin, st)
        binG = core.arrays['binG'].copy()
        st = np.random.default_rng(7)
        p = _call(core, T, x, want, ABI[formula], True, Cin[pat].copy(), st)
        assert np.array_equal(core.arrays['binG'], binG)
        if want & ABI['PCX_PASS_G']:
            assert np.array_equal(core.arrays['binGp'], binG[pat])
        made = formula != 'PCX_FORMULA_PASSES'
        packed_j = bool(want & ABI['PCX_WANT_PACKED'])
        for k, (dv, pv) in enumerate(zip(d[:3], p[:3])):
            if k == 0 and not want & ABI['PCX_PASS_C']:
                continue                                   # (C not made)
            if k == 1 and (not made or packed_j):
                continue                                   # (Geq not made)
            if k == 2 and (not made or not want & ABI['PCX_PASS_G']):
                continue                                   # (J not made)
            want_v = dv[:len(pat)] if (k == 2 and packed_j) else dv[pat]
            assert np.array_equal(pv, want_v, equal_nan=True), (k, formula, wkey)
        for dv, pv in zip(d[3:], p[3:]):
            assert np.array_equal(dv, pv, equal_nan=True)

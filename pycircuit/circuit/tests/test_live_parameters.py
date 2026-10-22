"""A write to a parameter that no per-circuit cache holds keeps the circuit's
set-up: the stamp plan, the PCNR view and its C record survive a source's
value written between two solves, and every analysis reads the new value.

The caches are keyed on one parameter epoch (`ParameterDict._epoch`); a
class names the parameters it only reads at evaluation time
(`_live_params`), and a write that changes only those -- or changes nothing
-- does not move it.  The class's `update` must be the declaring class's:
a subclass that overrides it is not live.
"""
import warnings

import numpy as np
import pytest

from pycircuit.circuit import circuit, gnd
from pycircuit.circuit import elements_hdl as eh
from pycircuit.circuit.analysis_ss import AC, Noise
from pycircuit.circuit.dcanalysis import DC
from pycircuit.circuit.elements import IS, VS, ISin, R, SubCircuit, VPulse, VSin
from pycircuit.circuit.transient import Transient
from pycircuit.utilities.param import ParameterDict

TK = circuit.numeric


def _cir(src_cls=VS, **src):
    """A source driving a diode through a resistor (the diode a PCNR
    participant whose kernels are bound, so a DC builds the PCNR view and
    its C record)."""
    c = SubCircuit(toolkit=TK)
    n1, n2 = c.add_nodes('1', '2')
    c['V'] = src_cls(n1, gnd, toolkit=TK, **src)
    c['R'] = R(n1, n2, r=1e3, toolkit=TK)
    c['D'] = eh.DiodeSpiceHdl(n2, gnd, IS=1e-14, toolkit=TK)
    return c


def _dc(c):
    with warnings.catch_warnings():
        warnings.simplefilter('ignore')
        return np.asarray(DC(c).solve().x, dtype=float)


def _kept(c):
    """The circuit's per-circuit set-up, by identity."""
    d = c.__dict__
    return (d.get('_stamp_plan'), d.get('_pcnr_view', (None, None))[1], d.get('_pcnr_c'))


@pytest.mark.parametrize('src_cls, src, name, value', [
    (VS, {'v': 1.0}, 'v', 2.0),
    (VSin, {'vo': 1.0, 'va': 0.0}, 'v', 2.0),
    (VPulse, {'v': 1.0}, 'v', 2.0),
], ids=['VS', 'VSin', 'VPulse'])
def test_a_source_value_written_between_two_dcs_keeps_the_set_up(src_cls, src, name, value):
    """The second DC runs on the first's plan, view and C record, and its
    answer is the circuit's built with the new value, to the bit."""
    c = _cir(src_cls, **src)
    _dc(c)
    kept = _kept(c)
    assert all(k is not None for k in kept), kept
    epoch = ParameterDict._epoch
    setattr(c['V'].ipar, name, value)
    assert ParameterDict._epoch is epoch
    got = _dc(c)
    assert all(a is b for a, b in zip(_kept(c), kept))
    assert got.tobytes() == _dc(_cir(src_cls, **{**src, name: value})).tobytes()


def test_a_written_source_value_reaches_the_transient_ac_and_noise():
    """The transient, AC and noise of a circuit whose source's `v`, `vac`,
    `phase` and `noisePSD` were written after a DC equal a fresh circuit's,
    to the bit."""
    new = {'v': 1.5, 'vac': 0.5, 'phase': 30.0, 'noisePSD': 1e-16}
    c = _cir(VS, v=1.0)
    _dc(c)
    for k, v in new.items():
        setattr(c['V'].ipar, k, v)
    ref = _cir(VS, **new)
    f = np.array([1e3, 1e6])

    def runs(cir):
        x = _dc(cir)
        tr = Transient(cir, toolkit=TK).solve(tend=1e-6, timestep=1e-8)
        ac = AC(cir, dcx=x).solve(f)
        n2 = cir.nodenames['2']
        nz = Noise(cir, inputsrc='V', outputnodes=(n2, gnd), dcx=x).solve(f)
        svn = nz['Svnout']
        return (np.asarray(tr.v('2', gnd).get_y(), dtype=float).tobytes(),
                np.asarray(ac.v('2', gnd).get_y()).tobytes(),
                np.asarray(getattr(svn, 'get_y', lambda: svn)()).tobytes())

    assert runs(c) == runs(ref)


def test_a_current_source_s_value_keeps_the_set_up():
    c = SubCircuit(toolkit=TK)
    n1 = c.add_node('1')
    c['I'] = IS(gnd, n1, i=1e-3, toolkit=TK)
    c['D'] = eh.DiodeSpiceHdl(n1, gnd, IS=1e-14, toolkit=TK)
    _dc(c)
    kept = _kept(c)
    c['I'].ipar.i = 2e-3
    got = _dc(c)
    assert all(a is b for a, b in zip(_kept(c), kept))
    ref = SubCircuit(toolkit=TK)
    m1 = ref.add_node('1')
    ref['I'] = IS(gnd, m1, i=2e-3, toolkit=TK)
    ref['D'] = eh.DiodeSpiceHdl(m1, gnd, IS=1e-14, toolkit=TK)
    assert got.tobytes() == _dc(ref).tobytes()


def test_a_written_waveform_parameter_keeps_the_set_up_and_reaches_the_transient():
    """A pulse's width written after a DC: the set-up kept, the transient
    the one of a pulse built with it (the waveform rebuilt by `update`)."""
    p = {'v1': 0.0, 'v2': 1.0, 'td': 1e-7, 'tr': 1e-8, 'tf': 1e-8, 'pw': 2e-7, 'per': 1e-6}
    c = _cir(VPulse, **p)
    _dc(c)
    kept = _kept(c)
    c['V'].ipar.pw = 4e-7
    _dc(c)
    assert all(a is b for a, b in zip(_kept(c), kept))

    def tran(cir):
        r = Transient(cir, toolkit=TK).solve(tend=1e-6, timestep=1e-8)
        return np.asarray(r.v('2', gnd).get_y(), dtype=float).tobytes()

    assert tran(c) == tran(_cir(VPulse, **{**p, 'pw': 4e-7}))


def test_a_baked_parameter_still_rebuilds():
    """A resistor's `r` is baked into the plan's constant stamps: its write
    moves the epoch and the answer is the new circuit's."""
    c = _cir(VS, v=1.0)
    _dc(c)
    plan = c.__dict__['_stamp_plan']
    epoch = ParameterDict._epoch
    c['R'].ipar.r = 2e3
    assert ParameterDict._epoch is not epoch
    got = _dc(c)
    assert c.__dict__['_stamp_plan'] is not plan
    ref = _cir(VS, v=1.0)
    ref['R'].ipar.r = 2e3
    assert got.tobytes() == _dc(ref).tobytes()


def test_a_participant_s_parameter_still_rebuilds_the_view():
    c = _cir(VS, v=1.0)
    _dc(c)
    view = c.__dict__['_pcnr_view'][1]
    c['D'].ipar.IS = 2e-14
    _dc(c)
    assert c.__dict__['_pcnr_view'][1] is not view


def test_a_subclass_that_overrides_update_is_not_live():
    """Its `update` may bake the value: a write moves the epoch."""
    class Baking(VS):
        def update(self, subject):
            super().update(subject)
            self._baked = self.iparv.v

    s = Baking(1, gnd, v=1.0, toolkit=TK)
    epoch = ParameterDict._epoch
    s.ipar.v = 2.0
    assert ParameterDict._epoch is not epoch
    assert s._baked == 2.0


def test_an_update_iparv_that_changes_nothing_moves_nothing():
    """`update_iparv` re-resolves every parameter; where no value (and no
    given flag) changes, the epoch holds and the set-up is kept."""
    c = _cir(VS, v=1.0)
    _dc(c)
    kept = _kept(c)
    epoch = ParameterDict._epoch
    c.update_iparv()
    assert ParameterDict._epoch is epoch
    _dc(c)
    assert all(a is b for a, b in zip(_kept(c), kept))


def test_a_given_flag_that_changes_alone_moves_the_epoch():
    """A parameter written to its default value changes what was GIVEN,
    which a compact model may branch on (`$param_given`): the epoch moves
    though no value does."""
    r = R(1, gnd, r=1e3, toolkit=TK)
    epoch = ParameterDict._epoch
    r.ipar.set(noisy=r.ipar.noisy)
    assert ParameterDict._epoch is not epoch


def test_a_copied_source_keeps_its_live_names():
    """`copy` of an element copies its parameter dicts; the copies keep the
    live names, so a copied circuit's sources keep the set-up too."""
    from copy import copy
    s = copy(VS(1, gnd, v=1.0, toolkit=TK))
    epoch = ParameterDict._epoch
    s.ipar.v = 2.0
    assert ParameterDict._epoch is epoch
    assert s.iparv.v == 2.0


def test_the_live_names_are_each_class_s_own():
    """Each class's live parameters are ones it declares; a waveform
    source's live set covers its waveform (built by `update`)."""
    for cls in (VS, IS, VSin, ISin, VPulse):
        names = {p.name for p in cls.instparams}
        assert set(cls._live_params) <= names, cls
    assert {'freq', 'phase', 'acphase'} <= set(VSin._live_params)

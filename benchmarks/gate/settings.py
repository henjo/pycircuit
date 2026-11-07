"""What a gate runs: the gold decks, the forced settings that take them through code they do not reach by
themselves, the switches that turn the C paths off, and the circuits at size.

A forcing sets a module's constants or a run's parameters; it moves no line of the simulator.  Each earns its place by
what it reached on the decks when measured on 866a0300 (roadmap 288): `mi6` (a Newton cut short: hand-backs, the fallback's limit,
the DC Newton's C) most, `bound` (the error test on the seed's factors whatever the chord's iterations, so no full test
comes before the first fallback: the C-made first numeric on 7 decks), `low` (every size gate at nothing, so packed sets,
pattern pages and the results' gathering run on small decks).  Forcing KLU alone, plain 'seed' and `low` with `bound`
added at most 2 lines and are not in the matrix.
"""

DECKS = ('4049osc', 'toronto', 'voter25', 'ab_integ', 'e1480', 'g1310', 'gm6', 'hussamp', 'mosrect', 'mux8', 'nand', 'ring',
         'slowlatch', 'rca', 'schmitecl', 'latch', 'opampal', 'gm1', 'gm2', 'gm3', 'gm17', 'mike2', 'rich3', 'todd3',
         'arom', 'gm19', 'jge')

#: (method, setting): every deck, the C paths against the Python's
MATRIX = (('gear2', 'default'), ('gear2', 'mi6'), ('gear2', 'bound'), ('gear2', 'low'),
          ('radau', 'default'), ('radau', 'mi6'), ('radau', 'low'),
          ('dc', 'default'))

#: the switches that take a run off its C stepping paths ('py'), and off every C path that has a Python ('off': slower)
PATHS = {
    'c': {},
    'py': {'gear2': ('PYCIRCUIT_LOOP_C', 'PYCIRCUIT_NEWTON_C'),
           'radau': ('PYCIRCUIT_TCLOOP_C', 'PYCIRCUIT_STAGELOOP_C'),
           'dc': ('PYCIRCUIT_PCNR_C', 'PYCIRCUIT_DCNEWTON_C')},
    'off': {'gear2': ('PYCIRCUIT_LOOP_C', 'PYCIRCUIT_NEWTON_C', 'PYCIRCUIT_LTE_C', 'PYCIRCUIT_CTRL_C', 'PYCIRCUIT_PRED_C',
                      'PYCIRCUIT_BRANCH_C', 'PYCIRCUIT_SOURCES_C'),
            'radau': ('PYCIRCUIT_TCLOOP_C', 'PYCIRCUIT_STAGELOOP_C', 'PYCIRCUIT_RADAU_C', 'PYCIRCUIT_RADAU_TC', 'PYCIRCUIT_LTE_C',
                      'PYCIRCUIT_CTRL_C', 'PYCIRCUIT_PRED_C', 'PYCIRCUIT_BRANCH_C', 'PYCIRCUIT_SOURCES_C'),
            'dc': ('PYCIRCUIT_PCNR_C', 'PYCIRCUIT_DCNEWTON_C')},
}

#: C against Python at size: each past the input's first edge, where the first fallback comes.  A shorter run compares
#: nothing the size gates or the fallback do (12 ns: `newton_c:first_numeric` 0).  53-173 s each on the Python's side.
SIZE_TWINS = (('chain', '4000', 'gear2', 1e-6), ('mesh', '63x63', 'gear2', 1e-7), ('mesh', '63x63', 'radau', 1e-7))

#: at size, the C's own paths further: reach and digests, tree against tree
SIZE_RUNS = (('chain', '4000', 'radau', 1e-6), ('achain', '4000', 'gear2', 1e-7), ('achain', '4000', 'radau', 1e-7),
             ('mesh', '45x45', 'radau', 1e-7), ('mesh', '90x90', 'radau', 3e-8))


def force_modules(setting):
    """A setting's module constants, before the run's objects are made."""
    parts = set(setting.split('+'))
    if 'bound' in parts:
        from pycircuit.circuit import _tran_newton, _tran_newton_c
        _tran_newton.LTE_SEED_ITERS = 1000
        _tran_newton_c._mods()['lte_seed'] = 1000
    if 'low' in parts:
        from pycircuit.circuit import (
            _cbuild,
            _tran_core,
            _tran_history,
            _tran_newton_c,
            _tran_radau,
            _tran_radau_c,
            _tran_tcloop_c,
        )
        from pycircuit.circuit.transient import Transient
        _tran_newton_c.SETS_FROM = 0
        _tran_tcloop_c.SETS_FROM = 0
        _tran_tcloop_c.PACKED = 'always'
        _tran_core.PATTERN_PAGES_FROM = 0
        _tran_radau.PATTERN_FROM = 0
        _tran_radau_c.OUT_BYTES = 1
        _tran_history.GATHER_MIN_BYTES = 0
        _tran_history.GATHER_TRIM_BYTES = 1 << 20
        _cbuild.INPLACE_MIN_ENTRIES = 0
        Transient._STRUCTURAL_KLU_MIN_N = property(lambda self: 8)


def force_run(setting, tr):
    """A setting's parameters on the transient `tr`."""
    parts = set(setting.split('+'))
    if parts & {'seed', 'bound'}:
        tr.par.lte_jacobian = 'seed'
    for p in parts:
        if p.startswith('mi') and p[2:].isdigit():
            tr.par.maxiter = int(p[2:])

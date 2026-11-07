"""THE GATE'S DECISION, apart from its measurement (`benchmarks/gate/verdict.py`, roadmap 288): a verdict has three outcomes, SAME, NOT REACHED and
DIFFERS, and each branch is tested here on rows made by hand, whatever the machine is doing.  "0 differing" with the code unreached read as a pass
twice on 2026-11-07: both shapes are pinned below.  And the coverage build (`_cbuild`'s `PYCIRCUIT_C_COVERAGE`): with it unset, every key is what it was."""
import importlib.util
import os
import subprocess
import sys

import numpy as np
import pytest

from pycircuit.circuit import _cbuild

ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), '..', '..', '..'))
## (loaded from its file, not imported: `benchmarks/` stays off `sys.path`)
_spec = importlib.util.spec_from_file_location('_gate_verdict', os.path.join(ROOT, 'benchmarks', 'gate', 'verdict.py'))
V = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(V)

GOOD = {'x': 'a', 'err': None, 'steps': 10, 'rej': 1, 'newton': 20, 'fb': 0, 'klu': {'factors': 1, 'refactors': 2}, 'counts': {'newton_c:first_numeric': 1}}


def _row(group, deck, path, **over):
    return {'group': group, 'deck': deck, 'path': path, 'row': dict(GOOD, **over)}


def test_held_memory_is_compared_where_both_sides_measured_it():
    """KLU's held bytes (`klumem`): a leak differs; a side that could not measure them (a tree without the reader, no KLU) is not a difference."""
    assert V.differing(dict(GOOD, kmem=100), dict(GOOD, kmem=140)) == {'kmem': (100, 140)}
    assert V.differing(dict(GOOD, kmem=100), dict(GOOD, kmem=None)) == {}
    assert V.differing(dict(GOOD, zmem=None), dict(GOOD, zmem=8)) == {}


def test_klu_s_own_memory_count_is_read_and_confirmed():
    """`klumem` reads `klu_common.memusage` at the offset the compiler gives on the installed header, trusted after KLU's own behaviour confirms it: a factor
    holds memory, its objects freed hold none."""
    from pycircuit._testing import klumem
    if klumem._offset() is None:
        pytest.skip('KLU, its header or a compiler absent')
    from pycircuit.circuit.linearsolver import KLUSolver
    sk = KLUSolver()
    assert klumem.memusage(sk) == 0
    sk.solve(np.array([[2.0, 1.0], [1.0, 3.0]]), np.ones(2), None)
    assert klumem.memusage(sk) > 0
    sk._free_numeric()
    sk._free_symbolic()
    assert klumem.memusage(sk) == 0
    assert klumem.memusage(object()) is None


def test_differing_names_what_differs_and_nothing_else():
    assert V.differing(GOOD, dict(GOOD)) == {}
    assert V.differing(GOOD, dict(GOOD, steps=11)) == {'steps': (10, 11)}
    assert V.differing(GOOD, dict(GOOD, klu={'factors': 0, 'refactors': 3})).keys() == {'klu'}
    assert V.differing(GOOD, dict(GOOD, counts={})) == {}, 'a counter is reach, not a result'


@pytest.mark.parametrize('a, b', [(None, GOOD), (GOOD, None), ({'failed': True}, GOOD), ({'failed': True}, {'failed': True})])
def test_a_missing_or_failed_side_differs_even_from_another_failure(a, b):
    assert V.differing(a, b).keys() == {'failed'}


def test_compare_pairs_by_group_and_deck():
    rows = [_row('g', 'a', 'c'), _row('g', 'a', 'py'), _row('g', 'b', 'c', x='b'), _row('g', 'b', 'py', x='c'), _row('h', 'a', 'c'), _row('h', 'a', 'py')]
    got = V.compare(rows)
    assert got['g']['runs'] == 2 and [d for d, _w in got['g']['differ']] == ['b']
    assert got['h'] == {'runs': 1, 'differ': []}


def test_a_deck_with_one_side_differs():
    assert [d for d, _w in V.compare([_row('g', 'a', 'c')])['g']['differ']] == ['a']


def test_a_group_with_neither_side_is_not_a_comparison():
    assert V.compare([_row('g', 'a', 'off')]) == {}


def test_with_side_keeps_only_what_can_be_compared_with_it():
    rows = [_row('g', 'a', 'c'), _row('g', 'a', 'off'), _row('g', 'b', 'c'), _row('h', 'a', 'c')]
    assert [(r['group'], r['deck'], r['path']) for r in V.with_side(rows, 'off')] == [('g', 'a', 'c'), ('g', 'a', 'off')]
    assert V.compare(V.with_side(rows, 'off'), right='off') == {'g': {'runs': 1, 'differ': []}}


def test_reach_sums_the_compared_path_and_shows_a_zero():
    rows = [_row('g', 'a', 'c'), _row('g', 'b', 'c'), _row('g', 'a', 'py', counts={'k': 100}), {'group': 'g', 'deck': 'c', 'path': 'c', 'row': {'failed': True}}]
    assert V.reach(rows, ['newton_c:first_numeric', 'k']) == {'newton_c:first_numeric': 2, 'k': 0}


@pytest.mark.parametrize('cls, sets', [('decks', ({1}, {1}, {1}, {1})), ('forced', (set(), {1}, {1}, {1})), ('size', (set(), set(), {1}, {1})),
                                       ('fixture', (set(), set(), set(), {1})), ('none', (set(), set(), set(), set()))])
def test_an_arm_takes_the_first_class_that_holds_it(cls, sets):
    assert V.arm_class(1, *sets) == cls


def test_changed_lines_reads_the_new_side():
    diff = """diff --git a/x/newton.c b/x/newton.c
--- a/x/newton.c
+++ b/x/newton.c
@@ -10,2 +12,3 @@ static long f(void)
+a
@@ -30 +33 @@
+b
@@ -40,3 +44,0 @@
diff --git a/x/gone.c b/x/gone.c
--- a/x/gone.c
+++ /dev/null
@@ -1,2 +0,0 @@
"""
    assert V.changed_lines(diff) == {'newton.c': {12, 13, 14, 33}}


def _changed(line_cls='decks', arm_cls='decks'):
    return V.changed_report({'n.c': {5, 6}}, {('n.c', 5)}, lambda k: line_cls, lambda k: [('n.c', k[1], 0)], lambda k: arm_cls)


def test_a_line_no_code_was_made_for_is_not_reported():
    got = _changed()
    assert got == {'lines': [(('n.c', 5), 'decks')], 'arms': [(('n.c', 5, 0), 'decks')]}


COMPARED = {'twin': {'runs': 27, 'differ': []}}


def test_same_needs_agreement_and_reach():
    assert V.verdict(COMPARED, {'k': 3}, _changed('forced', 'size')) == ('SAME', [])


def test_differs_outranks_everything_and_names_the_deck():
    out, why = V.verdict({'twin': {'runs': 27, 'differ': [('voter25', {'x': ('a', 'b')})]}}, {'k': 0}, _changed('none', 'none'))
    assert out == 'DIFFERS' and 'voter25' in why[0] and len(why) == 1


def test_nothing_compared_is_not_a_pass():
    assert V.verdict({})[0] == 'NOT REACHED'
    assert V.verdict({'twin': {'runs': 0, 'differ': []}})[0] == 'NOT REACHED'


def test_the_seed_twin_s_shape_reads_not_reached():
    """54 deck runs, 0 differing, the counter of the code under test at zero: 2026-11-07's first shape."""
    out, why = V.verdict({'twin': {'runs': 54, 'differ': []}}, {'newton_c:first_numeric': 0})
    assert out == 'NOT REACHED' and 'newton_c:first_numeric' in why[0]


def test_a_changed_line_only_a_fixture_takes_reads_not_reached():
    out, why = V.verdict(COMPARED, {}, _changed('fixture', 'decks'))
    assert out == 'NOT REACHED' and 'changed line n.c:5' in why[0] and 'fixture' in why[0]


def test_a_changed_arm_nothing_takes_reads_not_reached_on_a_reached_line():
    """A one-line `if (x) return BAIL;` that never fires: its line reached, its arm not."""
    out, why = V.verdict(COMPARED, {}, _changed('decks', 'none'))
    assert out == 'NOT REACHED' and why == ['changed arm n.c:5#0 not compared (none)']


def test_a_quiet_benchmark_is_named():
    rows = [{'name': 'chain', 'n': 4005, 'swing': 4}, {'name': 'achain', 'n': 4005, 'swing': 4000}, {'name': 'unmeasured', 'n': 10}]
    assert V.quiet(rows) == ['chain']


@pytest.mark.parametrize('ratios, op, bound, want', [([1.3, 1.2], '>=', 1.25, (1.25, True)), ([1.1, 1.2, 1.3], '>=', 1.25, (1.2, False)),
                                                     ([1.0], '<=', 1.1, (1.0, True)), ([], '>=', 1.0, (None, False))])
def test_a_timing_line_reads_the_median(ratios, op, bound, want):
    got = V.line_met(ratios, op, bound)
    assert got[1] == want[1] and (got[0] is None or abs(got[0] - want[0]) < 1e-12)


def test_what_has_no_twin_is_in_the_matrix_and_read_tree_against_tree():
    spec = importlib.util.spec_from_file_location('_gate_settings', os.path.join(ROOT, 'benchmarks', 'gate', 'settings.py'))
    settings = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(settings)
    assert settings.NO_TWIN <= set(settings.MATRIX)
    assert all(setting == 'default' or (method, setting) in settings.NO_TWIN or setting in ('mi6', 'bound', 'low') for method, setting in settings.MATRIX)


def test_the_coverage_build_is_off_and_moves_no_key_by_default(monkeypatch):
    """Unset, a runtime unit compiles with `CFLAGS` alone, so its key is the one a tree without the switch makes."""
    monkeypatch.setattr(_cbuild, 'COVERAGE_FLAGS', ())
    assert _cbuild._flags(False) == tuple(_cbuild.CFLAGS) == _cbuild._flags(True)
    key = _cbuild.source_key('int f(void);', False)
    monkeypatch.setattr(_cbuild, 'COVERAGE_FLAGS', ('--coverage',))
    assert _cbuild.source_key('int f(void);', False) != key


def test_the_coverage_build_takes_the_units_and_not_the_kernels(monkeypatch):
    monkeypatch.setattr(_cbuild, 'COVERAGE_FLAGS', ('--coverage',))
    assert '--coverage' in _cbuild._flags(False) and '--coverage' not in _cbuild._flags(True)
    kernel = _cbuild.source_key('int f(void);', True)
    monkeypatch.setattr(_cbuild, 'COVERAGE_FLAGS', ())
    assert _cbuild.source_key('int f(void);', True) == kernel


def test_the_coverage_build_writes_its_notes_and_counts_where_it_is_told(tmp_path):
    """A process on the switch: the unit's notes in its directory, its counts under `GCOV_PREFIX` at exit, with its own store."""
    notes = tmp_path / 'notes'
    notes.mkdir()
    env = dict(os.environ, PYCIRCUIT_C_COVERAGE=str(notes), PYCIRCUIT_HDL_CACHE_DIR=str(tmp_path / 'store'),
               GCOV_PREFIX=str(tmp_path / 'scope'), GCOV_PREFIX_STRIP=str(str(notes).count('/')))
    code = ('from pycircuit.circuit import _cbuild\n'
            "ffi, fn, key, cold, _s = _cbuild.load_unit('int pcx_t(int x) { return x > 0 ? x : -x; }\\n', 'int pcx_t(int x);', 'pcx_t')\n"
            'assert cold and fn(-3) == 3 and fn(4) == 4\n')
    subprocess.run([sys.executable, '-c', code], env=env, cwd=ROOT, check=True, timeout=120)
    assert any(p.suffix == '.gcno' for p in notes.iterdir())
    assert any(p.suffix == '.gcda' for p in (tmp_path / 'scope').iterdir())

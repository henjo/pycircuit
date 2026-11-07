"""A mutation sample's line changes (`benchmarks/gate/mutate.py`, roadmap 289), and the path counters a gate compares (`verdict.count_moves`).  Pure: on text and dictionaries."""
import importlib.util
import os
import re

import pytest

ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), '..', '..', '..'))


def _load(name):
    ## (loaded from its file, not imported: `benchmarks/` stays off `sys.path`)
    spec = importlib.util.spec_from_file_location('_gate_' + name, os.path.join(ROOT, 'benchmarks', 'gate', name + '.py'))
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


M = _load('mutate')
V = _load('verdict')


def _kinds(line):
    return {k: [x.strip() for kk, x in M.mutations(line) if kk == k] for k, _x in M.mutations(line)}


def test_a_comparison_flips_to_its_neighbour_and_a_pointer_is_left_alone():
    assert _kinds('    if (a < b) f(c);')['cmp'] == ['if (a <= b) f(c);']
    assert _kinds('    if (a >= b) f(c);')['cmp'] == ['if (a > b) f(c);']
    assert _kinds('    if (a == b) f(c);')['cmp'] == ['if (a != b) f(c);']
    assert 'cmp' not in _kinds('    s->x = t->y;')
    assert 'cmp' not in _kinds('    x = y >> 2;')


def test_a_condition_negates_only_where_it_closes_on_its_line():
    assert _kinds('    if (a && b) f(c);')['neg'] == ['if (!(a && b)) f(c);']
    assert 'neg' not in _kinds('    if (a && (b ||')


def test_the_other_kinds():
    assert _kinds('    if (a && b) f(c);')['logic'] == ['if (a || b) f(c);']
    assert _kinds('    x = y + z;')['arith'] == ['x = y - z;']
    assert _kinds('    x = -y;').get('arith') is None, 'a unary minus is not an operation to swap'
    assert _kinds('    x = y[2];')['const'] == ['x = y[3];']
    assert _kinds('    x = 1.5;')['const'] == ['x = 3.0;']
    assert _kinds('    return 0;')['ret'] == ['return 1;']


@pytest.mark.parametrize('line', ['    s->n++;', '    x = y;', '    foo(a, b);'])
def test_a_plain_statement_can_be_dropped(line):
    assert _kinds(line)['drop'] == [';']


@pytest.mark.parametrize('line', ['    double h;', '    return x;', '    long i = 0;', '    foo(a,', '    if (a) b();', '    } else {'])
def test_what_is_not_a_plain_statement_is_not_dropped(line):
    assert 'drop' not in _kinds(line)


@pytest.mark.parametrize('line', ['#define X (a < b)', '    /* a < b */', ' * a < b', '    // a < b', '    printf("a < b");', ''])
def test_a_line_with_nothing_to_read_gives_nothing(line):
    assert M.mutations(line) == []


def test_a_comment_on_the_line_is_not_changed():
    assert all('/* c < d */' in x for _k, x in M.mutations('    if (a < b) f(c);  /* c < d */'))
    assert len([k for k, _x in M.mutations('    x = y;  /* c < 3 */')]) == 1


def test_no_change_comes_twice():
    news = [x for _k, x in M.mutations('    return 0;')]
    assert len(news) == len(set(news))


def _source():
    return {'a.c': {i: f'    x{i} = y + {i};' for i in range(1, 301)}, 'b.c': {i: f'    if (p{i} < q) f();' for i in range(1, 101)}}


def test_a_draw_is_reproducible_and_spread_by_size():
    one, two = M.draw(_source(), 40, 7), M.draw(_source(), 40, 7)
    assert one == two and one != M.draw(_source(), 40, 8)
    assert one['a.c']['want'] == 30 and one['b.c']['want'] == 10


def test_a_small_unit_still_gets_its_least():
    assert M.draw({'a.c': {i: f'    x = {i};' for i in range(1, 1001)}, 'b.c': {1: '    x = 1;'}}, 10, 1)['b.c']['want'] == 3


def test_every_drawn_change_changes_its_line():
    for d in M.draw(_source(), 40, 3).values():
        assert all(m['old'] != m['new'] for m in d['drawn'])


def test_a_counter_that_moved_is_named_and_a_once_key_is_not():
    a, b = {'k': 1, 'once:x': 1, 'same': 2}, {'k': 2, 'once:x': 9, 'same': 2, 'new': 3}
    assert V.count_moves(a, b) == {'k': (1, 2), 'new': (0, 3)}


def test_a_stage_names_the_counters_it_moves_on_purpose():
    assert V.count_moves({'loop_c:bail:numeric': 2}, {}, re.compile('bail:numeric')) == {}
    assert V.count_moves({'loop_c:bail:other': 2}, {}, re.compile('bail:numeric')) == {'loop_c:bail:other': (2, 0)}


def test_counters_are_compared_where_both_ran():
    rows = [{'group': 'g', 'deck': 'a', 'path': 'new', 'row': {'counts': {'k': 1}}}, {'group': 'g', 'deck': 'a', 'path': 'base', 'row': {'counts': {'k': 2}}},
            {'group': 'g', 'deck': 'b', 'path': 'new', 'row': {'failed': True}}, {'group': 'g', 'deck': 'b', 'path': 'base', 'row': {'counts': {}}}]
    got = V.compare_counts(rows)
    assert got['g']['runs'] == 1 and got['g']['differ'] == [('a', {'k': (1, 2)})]

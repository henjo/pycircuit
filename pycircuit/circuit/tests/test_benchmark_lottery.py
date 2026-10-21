"""`benchmarks/lottery.py`'s reading (2026-10-21; doc/hdl_roadmap_260824.md,
162): a deck's own noise from a lottery that is alive -- an arm counts only
where its waveform differs from the base's, a deck with none is dead, and the
bound a line cites is the slowest live arm's time plus the slack."""
import os

import pytest


@pytest.fixture(autouse=True, scope='module')
def _benchmarks_on_path():
    """`benchmarks/` importable for this module's tests, and only for them."""
    with pytest.MonkeyPatch.context() as mp:
        mp.syspath_prepend(os.path.join(os.path.dirname(os.path.abspath(__file__)),
                                        os.pardir, os.pardir, os.pardir, 'benchmarks'))
        yield


def _run(deck, f, s, steps, x, rnd=1):
    return {'deck': deck, 'method': 'radau', 'f': f, 's': s, 'steps': steps, 'rej': 0, 'x': x,
            'round': rnd}


def test_an_arm_counts_only_where_it_is_alive():
    """Two rounds: the base and three arms, one of them dead (the base's
    digest); the best of the rounds per arm, the span over the base and the
    live arms only, the bound the slowest live arm's plus 3 %."""
    import lottery
    rows = []
    for rnd, jitter in ((1, 0.0), (2, 0.05)):
        rows += [_run('d', 1.0, 1.00 + jitter, 100, 'base', rnd),
                 _run('d', 0.99999, 1.10 + jitter, 104, 'a1', rnd),
                 _run('d', 1.00001, 0.90 + jitter, 97, 'a2', rnd),
                 _run('d', 1.00003, 5.00 + jitter, 100, 'base', rnd)]
    tab = lottery.table(rows)
    t = tab['d']
    assert (t['live'], t['arms'], t['dead'], t['rounds']) == (2, 3, False, 2)
    assert (t['base_s'], t['min_s'], t['max_s']) == (1.00, 0.90, 1.10)
    assert (t['base_steps'], t['min_steps'], t['max_steps']) == (100, 97, 104)
    assert lottery.bound(tab, 'd') == pytest.approx(1.10 * 1.03)


def test_a_deck_with_no_live_arm_is_dead_and_bound_by_its_base():
    import lottery
    rows = [_run('d', 1.0, 2.0, 50, 'same'), _run('d', 0.99999, 1.0, 50, 'same'),
            _run('d', 1.00001, 3.0, 50, 'same')]
    tab = lottery.table(rows)
    assert tab['d']['dead'] and tab['d']['live'] == 0
    assert lottery.bound(tab, 'd') == pytest.approx(2.0 * 1.03)


def test_a_run_that_is_not_deterministic_has_no_lottery_to_read():
    import lottery
    rows = [_run('d', 1.0, 1.0, 10, 'x', 1), _run('d', 1.0, 1.0, 10, 'y', 2)]
    with pytest.raises(ValueError, match='not deterministic'):
        lottery.table(rows)


def test_a_failed_run_is_left_out():
    import lottery
    rows = [_run('d', 1.0, 1.0, 10, 'b'), _run('d', 0.99999, 1.2, 11, 'a'),
            {'deck': 'd', 'method': 'radau', 'f': 1.00001, 'error': 'timeout'}]
    tab = lottery.table(rows)
    assert tab['d']['arms'] == 1 and tab['d']['live'] == 1

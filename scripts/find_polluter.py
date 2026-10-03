#!/usr/bin/env python
"""Find the test that pollutes another, from a worker's recorded order
(2026-10-03).

    python scripts/find_polluter.py 'path::victim' <gate> gw5
    python scripts/find_polluter.py 'path::victim' <gate>/replay/.pytest-replay-gw5.txt

Takes the tests the worker ran BEFORE the victim (`test_replay.read`) and
bisects them with detect-test-pollution, which runs pytest serially with the
candidate list and the victim until one test is left whose presence makes the
victim fail.  The victim must FAIL in the worker's order and pass alone: a
moved result that does not fail (a digest, a last bit) has to be made a
failure first -- the leak detector's invariant does that for a class bound to
C without its kernels.

The sub-runs are serial (`PYTEST_ADDOPTS=-n 0`: detect-test-pollution passes
no pytest arguments of its own), write no timing record, and run with the
leak detector OFF: in `fail` mode it would fail the polluter and restore the
state, and the victim would then pass.
"""
import os
import subprocess
import sys
import tempfile

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import test_replay


def main(argv=None):
    argv = list(sys.argv[1:] if argv is None else argv)
    if len(argv) < 2:
        print(__doc__)
        return 2
    victim, where = argv[0], argv[1]
    worker = argv[2] if len(argv) > 2 else None
    path = test_replay.replay_file(where, worker)
    order, _lines = test_replay.read(path, victim)
    print(f'{len(order) - 1} tests ran before {victim} in {path}', flush=True)
    ## (the victim is part of the list: detect-test-pollution takes the
    ## candidates as everything in it but the failing test)
    with tempfile.NamedTemporaryFile('w', suffix='.txt', delete=False) as tmp:
        tmp.write('\n'.join(order) + '\n')
    env = dict(os.environ, PYCIRCUIT_TEST_TIMINGS='0', PYCIRCUIT_LEAKS='off',
               PYTEST_ADDOPTS=(os.environ.get('PYTEST_ADDOPTS', '') + ' -n 0 -p no:cacheprovider').strip())
    exe = os.path.join(os.path.dirname(sys.executable), 'detect-test-pollution')
    try:
        return subprocess.call([exe, '--failing-test', victim, '--testids-file', tmp.name], env=env)
    finally:
        os.unlink(tmp.name)


if __name__ == '__main__':
    sys.exit(main())

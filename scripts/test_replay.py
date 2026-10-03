#!/usr/bin/env python
"""Replay one xdist worker's exact test order in ONE process (2026-10-03).

The recorded gate runs with `--replay-record-dir=<gate>/replay`, so every
worker leaves `.pytest-replay-gwN.txt`: JSON lines, each test twice (at its
start and at its finish), in the order that worker ran them.  An anomaly in
one worker of eight -- the G97/G98 case: a class bound to C with no kernels
-- is then reproduced by running that worker's order serially:

    python scripts/test_replay.py <gate>/replay/.pytest-replay-gw5.txt
    python scripts/test_replay.py <gate> gw5 --until 'path::test_name'
    python scripts/test_replay.py <gate> gw5 --list            # the node ids only

`--until NODEID` stops after that test (the victim).  Everything after the
options' `--` goes to pytest unchanged, so the gate's own plugins and
options can be given (`-- -p tran_recorder -p no:cacheprovider`).  The run
is serial (`-n 0`; pytest-replay takes one file serially), collects only the
files the replay names, and writes no timing record.
"""
import argparse
import json
import os
import subprocess
import sys
import tempfile


def replay_file(where, worker=None):
    """The replay file: `where` itself, or `<where>/replay/.pytest-replay-<worker>.txt`
    (or `<where>/.pytest-replay-<worker>.txt`)."""
    if worker is None:
        return where
    for cand in (os.path.join(where, 'replay', f'.pytest-replay-{worker}.txt'),
                 os.path.join(where, f'.pytest-replay-{worker}.txt')):
        if os.path.exists(cand):
            return cand
    raise SystemExit(f'no replay file for {worker} under {where}')


def read(path, until=None):
    """`(nodeids, lines)`: the node ids in the order the worker STARTED them,
    each once, cut after `until`; and the file's lines for those tests (both
    the start and the finish record, as pytest-replay reads them back)."""
    order, lines = [], []
    with open(path, encoding='utf-8') as f:
        recs = [json.loads(ln) for ln in f if ln.strip() and not ln.lstrip().startswith(('#', '//'))]
    for rec in recs:
        nid = rec['nodeid']
        if nid not in order:
            order.append(nid)
            if until is not None and nid == until:
                break
    if until is not None and until not in order:
        raise SystemExit(f'{until!r} is not in {path}')
    keep = set(order)
    lines = [json.dumps(r) for r in recs if r['nodeid'] in keep]
    return order, lines


def main(argv=None):
    argv = list(sys.argv[1:] if argv is None else argv)
    passthrough = []
    if '--' in argv:
        k = argv.index('--')
        argv, passthrough = argv[:k], argv[k + 1:]
    ap = argparse.ArgumentParser(description=__doc__.split('\n')[0])
    ap.add_argument('where', help='a replay file, or a gate directory with a worker name')
    ap.add_argument('worker', nargs='?', help='gwN')
    ap.add_argument('--until', help='stop after this node id (the victim)')
    ap.add_argument('--list', action='store_true', help='print the node ids and stop')
    a = ap.parse_args(argv)
    path = replay_file(a.where, a.worker)
    order, lines = read(path, a.until)
    if a.list:
        print('\n'.join(order))
        return 0
    files = sorted({nid.split('::', 1)[0] for nid in order})
    with tempfile.NamedTemporaryFile('w', suffix='.txt', delete=False) as tmp:
        tmp.write('\n'.join(lines) + '\n')
    env = dict(os.environ, PYCIRCUIT_TEST_TIMINGS='0')
    cmd = [sys.executable, '-m', 'pytest', *files, '--replay', tmp.name, '-n', '0', *passthrough]
    print(f'replaying {len(order)} tests from {path}' + (f' until {a.until}' if a.until else ''),
          flush=True)
    try:
        return subprocess.call(cmd, env=env)
    finally:
        os.unlink(tmp.name)


if __name__ == '__main__':
    sys.exit(main())

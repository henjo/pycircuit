# -*- coding: utf-8 -*-
# Copyright (c) 2026 Pycircuit Development Team
# See LICENSE for details.
"""Compare symbolic solve/representation backends on a fixed circuit set.

Developer tooling, deliberately *not* shipped inside the library (see
``doc/ddd_plan.md`` §3.6): the circuit fixtures live in
``pycircuit.circuit.benchmark_circuits`` because tests and the documentation
build both want them, while the timing/memory/isolation machinery here is only
of interest when running comparisons.

What it records per ``(system, backend)`` pair:

``size``
    Representation size in the backend's own unit -- ``count_ops`` for
    expression backends, vertex count for graph backends.  These units are *not*
    identical, which is why ``payload`` is reported beside them.
``payload``
    Total ``count_ops`` of the symbols/expressions each unit carries, so a
    backend with fewer, fatter units cannot look artificially compact.
``build_s`` / ``eval_s``
    Construction and per-sweep-point evaluation time, reported **separately**.
    Conflating them is how an earlier backend was adopted on a misleading
    figure: the "10x faster" symengine result was timing numeric evaluation
    rather than symbolic extraction.
``peak_rss_kb``
    Peak resident memory of the child process.
``status``
    ``ok`` / ``timeout`` / ``killed`` / ``error``.  Several backends *hang*
    rather than fail (GiNaC past dim ~16, SoE ``to_ratio`` for N>=5), so "did
    not finish within T" is a result, not an exception to swallow.

Every case runs in a subprocess under a wall-clock timeout and an address-space
cap, because benchmarks explore exactly the regime that once exhausted memory on
this machine.

Usage::

    python benchmarks/symbolic_bench.py                     # whole suite
    python benchmarks/symbolic_bench.py -s ladder -b soe    # filtered
    python benchmarks/symbolic_bench.py --validate          # harness self-check
"""

import argparse
import json
import os
import resource
import subprocess
import sys
import time

import numpy as np
import sympy

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), os.pardir))

from pycircuit.circuit import benchmark_circuits as bc


DEFAULT_TIMEOUT = 120.0
DEFAULT_MEM_MB = 4096

TEST_FREQ = 1e6


## ---------------------------------------------------------------- systems --

SYSTEM_FACTORIES = {
    'rc_ladder': bc.rc_ladder,
    'semi': bc.rc_ladder_semi_symbolic,
    'mfb': lambda: bc.mfb_filter(),
    'dense': bc.dense_symbolic_matrix,
    'legacy_soe': bc.legacy_soe_ladder,
}


def resolve(spec):
    """Turn ``'rc_ladder:8'`` into a `BenchSystem` (so children can rebuild it)."""
    if ':' in spec:
        name, arg = spec.split(':', 1)
        return SYSTEM_FACTORIES[name](int(arg))
    return SYSTEM_FACTORIES[spec]()


def suite_specs():
    """The fixed benchmark set, as specs."""
    specs = ['rc_ladder:%d' % n for n in bc.LADDER_SIZES]
    specs += ['semi:%d' % n for n in bc.LADDER_SIZES]
    specs += ['mfb']
    specs += ['dense:%d' % n for n in range(3, 7)]
    return specs


## --------------------------------------------------------------- backends --

def _numeric_reference(system):
    """Transfer function value at `TEST_FREQ`, by substituting then solving."""
    env = dict(system.params)
    env[system.s] = 1j * 2 * np.pi * TEST_FREQ
    n = system.A.rows
    An = np.array([[complex(system.A[i, j].subs(env)) for j in range(n)]
                   for i in range(n)], dtype=complex)
    bn = np.array([complex(system.b[i].subs(env)) for i in range(n)], dtype=complex)
    xn = np.linalg.solve(An, bn)
    return xn[system.out_index] / xn[system.in_index]


def _sub_env(system):
    env = dict(system.params)
    env[system.s] = 1j * 2 * np.pi * TEST_FREQ
    return env


def run_symbolic_poly(system):
    """Fraction-free ``N/D`` via :class:`SymbolicPolyToolkit`."""
    from pycircuit.circuit.toolkit import symbolic_poly
    t0 = time.perf_counter()
    num, den = symbolic_poly.linearsolver_num_den(system.A, system.b)
    build_s = time.perf_counter() - t0

    H = num[system.out_index] / num[system.in_index]
    size = int(sum(sympy.count_ops(n) for n in num) + sympy.count_ops(den))

    env = _sub_env(system)
    t0 = time.perf_counter()
    value = complex(H.subs(env))
    eval_s = time.perf_counter() - t0
    return dict(size=size, payload=size, build_s=build_s, eval_s=eval_s,
                value=value)


def run_soe(system):
    """Sequence-of-expressions: shared assignments, never inlined."""
    from pycircuit.circuit.soe import solve_soe
    t0 = time.perf_counter()
    sol = solve_soe(system.A, system.b)
    build_s = time.perf_counter() - t0

    assigns = sol.assignments
    x = sol.solution
    H = x[system.out_index] / x[system.in_index]
    size = int(sum(sympy.count_ops(e) for _, e in assigns) + sympy.count_ops(H))

    env = _sub_env(system)
    t0 = time.perf_counter()
    e = dict(env)
    for tsym, expr in assigns:
        e[tsym] = complex(expr.subs(e))
    value = complex(H.subs(e))
    eval_s = time.perf_counter() - t0
    return dict(size=size, payload=size, build_s=build_s, eval_s=eval_s,
                value=value)


def run_ginac(system):
    """GiNaC native ``N/D`` -- optional; absent unless the extension is built."""
    from pycircuit.circuit.toolkit import ginac_toolkit
    if ginac_toolkit is None:
        raise RuntimeError('ginac_toolkit unavailable')
    t0 = time.perf_counter()
    num, den = ginac_toolkit.linearsolver_num_den(system.A, system.b)
    build_s = time.perf_counter() - t0

    H = num[system.out_index] / num[system.in_index]
    size = int(sum(sympy.count_ops(n) for n in num) + sympy.count_ops(den))

    env = _sub_env(system)
    t0 = time.perf_counter()
    value = complex(H.subs(env))
    eval_s = time.perf_counter() - t0
    return dict(size=size, payload=size, build_s=build_s, eval_s=eval_s,
                value=value)


def run_ddd(system):
    """Determinant decision diagram: shared graph, never expanded."""
    from pycircuit.circuit.ddd import ddd_cramer
    import sympy as _sp

    out, inp = system.out_index, system.in_index
    t0 = time.perf_counter()
    den, nums = ddd_cramer(system.A, system.b, indices=[out, inp])
    build_s = time.perf_counter() - t0

    ## The transfer function is a ratio of two numerator diagrams -- the shared
    ## determinant cancels -- so the denominator is not part of the cost of an
    ## x_out/x_in query.  Report what the query actually needs.
    graphs = [nums[out]] if out == inp else [nums[out], nums[inp]]
    size = sum(g.size for g in graphs)
    payload = int(sum(_sp.count_ops(g.matrix[i, j])
                      for g in graphs for (i, j) in g.entries()))

    env = _sub_env(system)
    t0 = time.perf_counter()
    value = complex(nums[out].eval(env)) / complex(nums[inp].eval(env))
    eval_s = time.perf_counter() - t0
    return dict(size=size, payload=payload, build_s=build_s, eval_s=eval_s,
                value=value, den_size=den.size)


BACKENDS = {
    'symbolic_poly': run_symbolic_poly,
    'soe': run_soe,
    'ginac': run_ginac,
    'ddd': run_ddd,
}


## ------------------------------------------------------------ measurement --

def measure(spec, backend):
    """Run one case **in process**.  Returns a JSON-serialisable dict."""
    system = resolve(spec)
    out = dict(system=spec, backend=backend, dim=system.dim)
    try:
        res = BACKENDS[backend](system)
    except Exception as exc:                     # noqa: BLE001 - reported, not raised
        out.update(status='error', error='%s: %s' % (type(exc).__name__, exc))
        return out

    value = res.pop('value')
    try:
        ref = _numeric_reference(system)
        err = abs(value - ref) / max(1.0, abs(ref))
        out['correct'] = bool(err < 1e-9)
        out['rel_err'] = float(err)
    except Exception as exc:                     # noqa: BLE001
        out['correct'] = None
        out['rel_err'] = None
        out['error'] = 'reference failed: %s' % exc

    out.update(res)
    out['status'] = 'ok'
    out['peak_rss_kb'] = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
    return out


def run_isolated(spec, backend, timeout=DEFAULT_TIMEOUT, mem_mb=DEFAULT_MEM_MB):
    """Run one case in a subprocess under a timeout and an address-space cap.

    A backend that exhausts memory takes down only the child, and one that hangs
    is recorded as ``timeout`` rather than stalling the suite.
    """
    cmd = [sys.executable, os.path.abspath(__file__),
           '--child', '--system', spec, '--backend', backend]
    env = dict(os.environ)
    env.setdefault('MPLBACKEND', 'Agg')
    repo = os.path.join(os.path.dirname(os.path.abspath(__file__)), os.pardir)
    env['PYTHONPATH'] = os.path.abspath(repo) + os.pathsep + env.get('PYTHONPATH', '')

    def limit():
        cap = mem_mb * 1024 * 1024
        resource.setrlimit(resource.RLIMIT_AS, (cap, cap))

    t0 = time.perf_counter()
    try:
        proc = subprocess.run(cmd, capture_output=True, text=True, env=env,
                              timeout=timeout, preexec_fn=limit)
    except subprocess.TimeoutExpired:
        return dict(system=spec, backend=backend, status='timeout',
                    wall_s=time.perf_counter() - t0)
    wall = time.perf_counter() - t0

    for line in reversed(proc.stdout.splitlines()):
        if line.startswith('{'):
            rec = json.loads(line)
            rec['wall_s'] = wall
            return rec
    status = 'killed' if proc.returncode != 0 else 'error'
    return dict(system=spec, backend=backend, status=status,
                returncode=proc.returncode, wall_s=wall,
                error=(proc.stderr or '').strip()[-400:])


## -------------------------------------------------------------- reporting --

_COLS = ('system', 'backend', 'dim', 'status', 'size', 'payload',
         'build_s', 'eval_s', 'peak_rss_kb', 'correct')


def format_table(rows):
    def cell(r, c):
        v = r.get(c)
        if v is None:
            return '-'
        if isinstance(v, float):
            return '%.4g' % v
        return str(v)

    widths = {c: max(len(c), *(len(cell(r, c)) for r in rows)) for c in _COLS} \
        if rows else {c: len(c) for c in _COLS}
    out = ['  '.join(c.ljust(widths[c]) for c in _COLS),
           '  '.join('-' * widths[c] for c in _COLS)]
    for r in rows:
        out.append('  '.join(cell(r, c).ljust(widths[c]) for c in _COLS))
    return '\n'.join(out)


def validate():
    """Harness self-check: reproduce the published SoE op counts.

    The SoE table in ``doc/src/circuit/soe_symbolic.rst`` reports 73/157/241/325
    operations at N=4/8/12/16 on the hand-built ladder.  If this harness cannot
    reproduce our own published numbers it cannot referee anyone else's, so this
    gates Stage B.
    """
    expected = {4: 73, 8: 157, 12: 241, 16: 325}
    ok = True
    print('Harness validation -- SoE op counts on the legacy hand-built ladder')
    print('  N   measured  expected  match')
    for N in sorted(expected):
        rec = measure('legacy_soe:%d' % N, 'soe')
        got = rec.get('size')
        match = (got == expected[N])
        ok = ok and match and rec.get('correct') is not False
        print('  %-3d %-9s %-9d %s' % (N, got, expected[N], match))
    print('\nVALIDATION %s' % ('PASSED' if ok else 'FAILED'))
    return 0 if ok else 1


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument('-s', '--system', action='append',
                    help='system spec (repeatable), e.g. rc_ladder:8')
    ap.add_argument('-b', '--backend', action='append', choices=sorted(BACKENDS),
                    help='backend (repeatable); default all')
    ap.add_argument('-t', '--timeout', type=float, default=DEFAULT_TIMEOUT)
    ap.add_argument('-m', '--mem-mb', type=int, default=DEFAULT_MEM_MB)
    ap.add_argument('--json', help='write raw records here')
    ap.add_argument('--validate', action='store_true',
                    help='self-check against published SoE numbers, then exit')
    ap.add_argument('--child', action='store_true', help=argparse.SUPPRESS)
    args = ap.parse_args(argv)

    if args.child:
        print(json.dumps(measure(args.system[0], args.backend[0]), default=str))
        return 0

    if args.validate:
        return validate()

    specs = args.system or suite_specs()
    backends = args.backend or sorted(BACKENDS)

    rows = []
    for spec in specs:
        for backend in backends:
            rec = run_isolated(spec, backend, args.timeout, args.mem_mb)
            rows.append(rec)
            print('  ran %-22s %-14s %s' % (spec, backend, rec.get('status')),
                  file=sys.stderr)

    print()
    print(format_table(rows))
    if args.json:
        with open(args.json, 'w') as fh:
            json.dump(rows, fh, indent=2, default=str)
    return 0


if __name__ == '__main__':
    sys.exit(main())

"""Large circuits: import them and measure where pycircuit's dense path
ends (the SPICE benchmark plan's stage 5; the sparse path is its own plan).

    python benchmarks/large_circuits.py ibm [--fraction K] [--merge-shorts]
    python benchmarks/large_circuits.py mesh [--max-n N] [--steps S]
    python benchmarks/large_circuits.py footprint N

`ibm`: the IBM power grid ibmpg1t (fetched, `benchmarks/fetch_spice_suite.py`)
read and imported -- the first 1/K of its element lines with `--fraction K`
(16, 4: a part, to see how the cost grows before the whole) -- its
parse and build times, peak RSS, element counts, unknowns;
`--merge-shorts` unifies the nodes 0 V sources join (the importer's
`merge_shorts`).  No transient is run: the record says whether the guard
would refuse one.

`mesh`: a synthetic power grid in the IBM files' form (`rc_mesh`: grid
resistors, a capacitor to ground at every node, R-L pads to a 1.8 V supply,
PULSE current sinks), written as SPICE and imported, on a ladder of sizes
-- 20 fixed steps each, in a process of its own: build and step times,
the step's split by profile (assembly, the linear solve and its
conversion, sources, the rest), the solver chosen, peak RSS -- up to the
size the guard refuses.

The guard: a step holds dense n x n arrays (the assembled matrices, the
Jacobian, the dense-to-sparse conversion's work); their bytes per n^2 are
fitted from the slope of the ladder's last two rungs' RSS, x 1.25
(`footprint`), and any run estimated past half of MemAvailable is refused
rather than started.  ⚠ The guard is an estimate, not the protection: each
measurement runs in a cgroup capped at half of MemAvailable with no swap
(`_child`), since an earlier version, guard only, took the shared box down
(2026-10-07).
"""
import bz2
import cProfile
import json
import os
import pstats
import resource
import subprocess
import sys
import time
import warnings

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from pycircuit._testing import benchdata

#: bytes per unknown squared a step's dense arrays hold: the mesh ladder's
#: slope between its two largest rungs (n = 4681 -> 9226, 1736 -> 7201 MB
#: peak RSS, 2026-10-08; 69 between the two below it -- it grows); the
#: ladder refits it as it climbs
BYTES_PER_N2 = 91


def mem_available():
    """MemAvailable, bytes (Linux)."""
    with open('/proc/meminfo') as fh:
        for line in fh:
            if line.startswith('MemAvailable:'):
                return int(line.split()[1]) * 1024
    raise RuntimeError('no MemAvailable in /proc/meminfo')


def footprint(n, per_n2=BYTES_PER_N2):
    """A step's dense arrays at `n` unknowns, bytes."""
    return per_n2 * n * n


def guard(n, per_n2=BYTES_PER_N2):
    """None where a run at `n` unknowns fits in half of MemAvailable, else
    why not."""
    need, have = footprint(n, per_n2), mem_available()
    if need > have / 2:
        return (f'{n} unknowns: a step holds ~{need / 1e9:.1f} GB dense, past half of the '
                f'{have / 1e9:.1f} GB available')
    return None


def rc_mesh(nx, ny, pad_every=8, sink_every=3):
    """A power grid of `nx` x `ny` nodes in the IBM files' form, as SPICE
    text: 0.25 ohm grid resistors, 10 fF to ground at each node, a pad
    (0.25 ohm, 1 nH, a 1.8 V source) every `pad_every` nodes in each
    direction, a PULSE current sink every `sink_every`."""
    out = [f'* rc_mesh {nx} x {ny}']
    for i in range(nx):
        for j in range(ny):
            n = f'n_{i}_{j}'
            if i + 1 < nx:
                out.append(f'rh_{i}_{j} {n} n_{i + 1}_{j} 0.25')
            if j + 1 < ny:
                out.append(f'rv_{i}_{j} {n} n_{i}_{j + 1} 0.25')
            out.append(f'c_{i}_{j} {n} 0 10f')
            if i % pad_every == 0 and j % pad_every == 0:
                out += [f'rp_{i}_{j} {n} x_{i}_{j} 0.25', f'lp_{i}_{j} x_{i}_{j} y_{i}_{j} 1n',
                        f'vp_{i}_{j} y_{i}_{j} 0 1.8']
            if (i * ny + j) % sink_every == 0:
                td = 1e-11 * ((i * 7 + j * 13) % 10)
                out.append(f'is_{i}_{j} {n} 0 2e-5 pulse(2e-05, 5e-3, {td:g},  1e-10,  1e-10,  '
                           '1e-11,  3e-09)')
    out += ['.tran 1e-11 2e-10', '.end']
    return '\n'.join(out) + '\n'


def _rss_mb():
    return resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 1024


def _import(path, merge_shorts=False):
    from pycircuit.circuit import circuit, spice_import
    circuit.default_toolkit = circuit.numeric
    return spice_import.import_netlist(path, strict=False, **(
        {'merge_shorts': True} if merge_shorts else {}))


#: the step's pieces: a profile's functions, by the module path they live in
_PIECES = {'assembly': ('_stamp_plan', 'circuit.py', '_tran_core', 'analysis.py'),
           'linear solve': ('linearsolver', 'scipy', 'numpy/linalg'),
           'sources': ('func.py', 'elements.py')}


def _split(prof):
    """The profile's self time by piece (`_PIECES`; the rest 'other')."""
    st = pstats.Stats(prof)
    out = dict.fromkeys(list(_PIECES) + ['other'], 0.0)
    for (fname, _line, _name), (_cc, _nc, tt, _ct, _callers) in st.stats.items():
        piece = next((p for p, keys in _PIECES.items() if any(k in fname for k in keys)),
                     'other')
        out[piece] += tt
    return out


def mesh_one(nx, ny, steps):
    """One mesh: its record."""
    import numpy as np

    from pycircuit.circuit import _paths
    path = os.path.join(benchdata.cache_root(), 'runs', 'mesh', f'mesh_{nx}x{ny}.cir')
    os.makedirs(os.path.dirname(path), exist_ok=True)
    with open(path, 'w') as fh:
        fh.write(rc_mesh(nx, ny))
    t0 = time.perf_counter()
    imp = _import(path)
    rec = {'nx': nx, 'ny': ny, 'build_s': time.perf_counter() - t0,
           'elements': len(imp.elements), 'unknowns': imp.circuit.n}
    rec['per_element_us'] = rec['build_s'] / rec['elements'] * 1e6
    tr, kw = imp.transient()
    # the branch check timed on its own (its structural rank, once a run;
    # one per-step screen), and ON in the solve, as a user runs it: its dense
    # SVDs were 93 % of the n = 4681 run until 5608aab4 and read as 'linear
    # solve'; c2b8f948 keeps a constant C's screen verdict
    x0 = np.zeros(imp.circuit.n)
    t0 = time.perf_counter()
    rec['branch_rank'] = tr._branch_structural_rank()[0]
    rec['branch_rank_s'] = time.perf_counter() - t0
    t0 = time.perf_counter()
    tr._branch_screen(x0)
    rec['branch_screen_ms'] = (time.perf_counter() - t0) * 1e3
    before = _paths.snapshot()
    prof = cProfile.Profile()
    t0 = time.perf_counter()
    with warnings.catch_warnings():
        warnings.simplefilter('ignore')
        prof.enable()
        res = tr.solve(tend=steps * kw['timestep'], timestep=kw['timestep'],
                       fixed_timestep=True)
        prof.disable()
    rec['solve_s'] = time.perf_counter() - t0
    rec['step_ms'] = rec['solve_s'] / steps * 1e3
    rec['newton'] = tr.statistics.newton_iterations
    rec['split'] = _split(prof)
    rec['dc_op_s'] = next((ct for (_f, _l, name), (_cc, _nc, _tt, ct, _c)
                           in pstats.Stats(prof).stats.items()
                           if name == '_solve_operating_point'), None)
    rec['solver'] = sorted(k for k in _paths.since(before) if k.startswith('solver.auto:'))
    x = np.asarray(res.x)
    rec['finite'] = bool(np.isfinite(x).all())
    rec['rss_mb'] = _rss_mb()
    return rec


def ibm(fraction, merge_shorts):
    """ibmpg1t (or its first 1/`fraction`): read, import, count."""
    import numpy as np

    from pycircuit.utilities import spicenetlist
    src = benchdata.spice_data('ibmpg1t.spice.bz2', 'ibm_pg')
    if src is None:
        sys.exit('ibmpg1t is not fetched (benchmarks/fetch_spice_suite.py)')
    with bz2.open(src, 'rt') as fh:
        lines = fh.read().splitlines()
    head = [ln for ln in lines if ln.startswith(('*', '.'))]
    body = [ln for ln in lines if ln and not ln.startswith(('*', '.'))]
    body = body[:len(body) // fraction]
    path = os.path.join(benchdata.cache_root(), 'runs', 'ibm', f'ibmpg1t_1of{fraction}.spice')
    os.makedirs(os.path.dirname(path), exist_ok=True)
    with open(path, 'w') as fh:
        fh.write('* ibmpg1t' + (f', its first 1/{fraction}' if fraction > 1 else '') + '\n')
        fh.write('\n'.join(body + [h for h in head if not h.startswith('*')]) + '\n')
    rec = {'fraction': fraction, 'merge_shorts': merge_shorts, 'lines': len(body)}
    t0 = time.perf_counter()
    net = spicenetlist.read(path)
    rec['read_s'] = time.perf_counter() - t0
    kinds = {}
    for c in net.cards:
        kinds[c.kind] = kinds.get(c.kind, 0) + 1
    rec['kinds'] = kinds
    t0 = time.perf_counter()
    imp = _import(path, merge_shorts)
    rec['import_s'] = time.perf_counter() - t0
    rec['per_element_us'] = rec['import_s'] / max(len(imp.elements), 1) * 1e6
    rec['elements'] = len(imp.elements)
    n = rec['unknowns'] = imp.circuit.n
    rec['nodes'] = len(imp.circuit.nodes)
    rec['rss_mb'] = _rss_mb()
    rec['refused'] = guard(n)
    rec['footprint_gb'] = footprint(n) / 1e9
    del np
    return rec


def _child(argv):
    """Run one measurement in a fresh process; its JSON record.  The child
    runs in a cgroup capped at half of MemAvailable with no swap, so an
    estimate that is wrong kills the child, not the box (the guard alone let
    a ladder rung take the machine down, 2026-10-07)."""
    cap = mem_available() // 2
    cmd = [sys.executable, os.path.abspath(__file__)] + argv
    cmd = ['systemd-run', '--user', '--scope', '-q', '-p', f'MemoryMax={cap}',
           '-p', 'MemorySwapMax=0'] + cmd
    p = subprocess.run(cmd, capture_output=True, text=True, check=False)
    line = p.stdout.strip().splitlines()[-1] if p.stdout.strip() else ''
    try:
        return json.loads(line)
    except ValueError:
        if p.returncode in (-9, 137):
            return {'error': f'killed at the {cap / 1e9:.1f} GB cap'}
        return {'error': (p.stderr.strip().splitlines() or ['?'])[-1][:300]}


def main(argv):
    if not argv:
        sys.exit(__doc__)
    if argv[0] == '--mesh-one':
        nx, ny, steps = (int(a) for a in argv[1:4])
        print(json.dumps(mesh_one(nx, ny, steps)))
        return 0
    if argv[0] == '--ibm-one':
        print(json.dumps(ibm(int(argv[1]), argv[2] == '1')))
        return 0
    if argv[0] == 'footprint':
        n = int(argv[1])
        print(f'{n} unknowns: ~{footprint(n) / 1e9:.2f} GB dense; '
              f'{guard(n) or "within half of MemAvailable"}')
        return 0
    if argv[0] == 'ibm':
        fraction = int(argv[argv.index('--fraction') + 1]) if '--fraction' in argv else 1
        rec = _child(['--ibm-one', str(fraction), '1' if '--merge-shorts' in argv else '0'])
        print(json.dumps(rec, indent=1))
        return 0
    if argv[0] == 'mesh':
        max_n = int(argv[argv.index('--max-n') + 1]) if '--max-n' in argv else 1 << 20
        steps = int(argv[argv.index('--steps') + 1]) if '--steps' in argv else 20
        side = 16
        per_n2, last = BYTES_PER_N2, None
        while True:
            n_est = side * side + 3 * (((side + 7) // 8) ** 2)
            if n_est > max_n:
                break
            why = guard(n_est, per_n2)
            if why:
                print(f'refused: {why}', flush=True)
                break
            rec = _child(['--mesh-one', str(side), str(side), str(steps)])
            print(json.dumps(rec), flush=True)
            if 'error' in rec:
                break
            # the guard's bytes per n^2 from the SLOPE between the last two
            # rungs, x 1.25: the max over rungs over the first one's RSS read
            # the small rungs' interpreter noise (135 B/n^2 where the slope is
            # 69) and refused n = 9081, which peaks at 5.4 GB (2026-10-08)
            if last is not None and rec['unknowns'] > last['unknowns']:
                fit = ((rec['rss_mb'] - last['rss_mb']) * 2**20
                       / (rec['unknowns'] ** 2 - last['unknowns'] ** 2))
                per_n2 = max(BYTES_PER_N2, 1.25 * fit)
            last = rec
            side = int(side * 1.41421356 + 0.5)
        return 0
    sys.exit(__doc__)


if __name__ == '__main__':
    sys.exit(main(sys.argv[1:]))

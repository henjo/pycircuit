"""The benchmark discipline every timing script here shares (2026-10-03;
robust timing, stage 5 of the testing plan).

A measurement on this box is only worth what its conditions were: it is a
hybrid i7-12800HX (P-cores at 4.7-4.8 GHz, E-cores at 3.4 GHz), the CPU
governor is `powersave` with no isolated cores, and several sessions share
it -- their test suites, their benchmarks.  So:

* `lock()` -- one benchmark at a time across every session that uses this
  module, and never beside a test suite of this repository: benchmarks take
  `~/.cache/pycircuit/bench.lock` EXCLUSIVELY, the root conftest takes it
  SHARED for a whole test session.  Waiting is printed with the holder.
* `pin_threads()` -- the BLAS/OpenMP thread variables to 1, which has to
  happen BEFORE numpy is imported (a pool opened before it keeps its 24
  threads); `step_machinery`'s `pss` case ran on 24 threads until it did.
* `pin_cpu()` -- the process on one P-core (CPU 8 unless
  `PYCIRCUIT_BENCH_CPU`), so parent and child never land on different core
  types or migrate.
* `Idle` -- the share of time CPU 8's hyperthread sibling (CPU 9) and the
  other P-cores were busy during a timed interval, from `/proc/stat`; a
  round measured while the sibling was busy is discarded and re-run.
* `stamp()` -- the conditions, recorded with every result: commit, load
  average, governor, CPU.

What it cannot do without root (documented, never required): set the
governor to `performance` (`cpupower frequency-set -g performance`), isolate
cores.  Interleaving parent and child per round and reading PAIRED ratios is
what makes the result robust to the drift those would remove.
"""
import contextlib
import os
import subprocess
import sys
import time

LOCK_PATH = os.path.join(os.path.expanduser('~'), '.cache', 'pycircuit', 'bench.lock')
CPU = int(os.environ.get('PYCIRCUIT_BENCH_CPU', '8'))
THREAD_VARS = ('OMP_NUM_THREADS', 'OPENBLAS_NUM_THREADS', 'MKL_NUM_THREADS',
               'BLIS_NUM_THREADS', 'NUMEXPR_NUM_THREADS')


def pin_threads():
    """One BLAS/OpenMP thread.  Must run before numpy is imported (a
    second call, after an earlier one pinned them, is a no-op)."""
    if all(os.environ.get(v) == '1' for v in THREAD_VARS):
        return
    if 'numpy' in sys.modules:
        raise RuntimeError('pin_threads() after numpy was imported: its BLAS pool is open')
    for v in THREAD_VARS:
        os.environ[v] = '1'


def no_aslr():
    """The command prefix that runs a child without address-space layout
    randomisation (`setarch <arch> -R`, no root needed), or [] where it is
    not available."""
    import platform
    import shutil
    if shutil.which('setarch') is None:
        return []
    return ['setarch', platform.machine(), '-R']


def pin_cpu(cpu=None):
    """This process (and its future threads) on one CPU."""
    cpu = CPU if cpu is None else cpu
    try:
        os.sched_setaffinity(0, {cpu})
    except (AttributeError, OSError) as e:
        print(f'_bench: could not pin to CPU {cpu}: {e}', file=sys.stderr)


def sibling(cpu=None):
    cpu = CPU if cpu is None else cpu
    try:
        with open(f'/sys/devices/system/cpu/cpu{cpu}/topology/thread_siblings_list') as f:
            txt = f.read().strip()
    except OSError:
        return None
    ids = set()
    for part in txt.split(','):
        a, _, b = part.partition('-')
        ids.update(range(int(a), int(b or a) + 1))
    ids.discard(cpu)
    return min(ids) if ids else None


def _cpu_times():
    out = {}
    with open('/proc/stat') as f:
        for ln in f:
            if ln.startswith('cpu') and ln[3].isdigit():
                p = ln.split()
                vals = [int(v) for v in p[1:]]
                idle = vals[3] + vals[4]                 # idle + iowait
                out[int(p[0][3:])] = (sum(vals), idle)
    return out


class Idle:
    """Busy shares of the CPUs over an interval: `with Idle() as w: ...;
    w.busy[9]` -- or `w.sibling_busy` for the timed CPU's sibling."""

    def __init__(self, cpu=None):
        self.cpu = CPU if cpu is None else cpu
        self.sib = sibling(self.cpu)
        self.busy = {}

    def __enter__(self):
        self.t0 = _cpu_times()
        return self

    def __exit__(self, *exc):
        t1 = _cpu_times()
        for c, (tot1, idle1) in t1.items():
            tot0, idle0 = self.t0.get(c, (tot1, idle1))
            d = tot1 - tot0
            self.busy[c] = 0.0 if d <= 0 else 1.0 - (idle1 - idle0) / d
        return False

    @property
    def sibling_busy(self):
        return self.busy.get(self.sib, 0.0) if self.sib is not None else 0.0


@contextlib.contextmanager
def lock(exclusive=True, what=''):
    """The benchmark lock (see the module note); prints who holds it while
    waiting.  Yields nothing; a platform without `fcntl` runs unlocked."""
    try:
        import fcntl
    except ImportError:                                  # pragma: no cover
        yield
        return
    os.makedirs(os.path.dirname(LOCK_PATH), exist_ok=True)
    fd = os.open(LOCK_PATH, os.O_RDWR | os.O_CREAT, 0o644)
    mode = fcntl.LOCK_EX if exclusive else fcntl.LOCK_SH
    try:
        try:
            fcntl.flock(fd, mode | fcntl.LOCK_NB)
        except OSError:
            holder = _holders()
            print(f'_bench: waiting for the benchmark lock ({holder})', file=sys.stderr, flush=True)
            t0 = time.time()
            fcntl.flock(fd, mode)
            print(f'_bench: got the lock after {time.time() - t0:.0f} s', file=sys.stderr, flush=True)
        _note_holder(exclusive, what)
        try:
            yield
        finally:
            _drop_holder()
            fcntl.flock(fd, fcntl.LOCK_UN)
    finally:
        os.close(fd)


def _holder_file():
    return os.path.join(os.path.dirname(LOCK_PATH), 'bench.lock.holders',
                        f'{os.getpid()}')


def _note_holder(exclusive, what):
    try:
        os.makedirs(os.path.dirname(_holder_file()), exist_ok=True)
        with open(_holder_file(), 'w') as f:
            f.write(f"{'benchmark' if exclusive else 'test suite'} pid {os.getpid()} "
                    f"{what or ' '.join(sys.argv)[:120]}\n")
    except OSError:
        pass


def _drop_holder():
    try:
        os.unlink(_holder_file())
    except OSError:
        pass


def _holders():
    d = os.path.dirname(_holder_file())
    out = []
    try:
        for name in os.listdir(d):
            try:
                os.kill(int(name), 0)
            except (OSError, ValueError):
                continue                                 # a stale note
            with open(os.path.join(d, name)) as f:
                out.append(f.read().strip())
    except OSError:
        pass
    return '; '.join(out) or 'holder unknown'


def stamp():
    """The conditions of a measurement."""
    root = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
    try:
        commit = subprocess.check_output(['git', 'rev-parse', '--short', 'HEAD'], cwd=root,
                                         stderr=subprocess.DEVNULL, text=True).strip()
    except (OSError, subprocess.CalledProcessError):
        commit = 'unknown'
    try:
        with open(f'/sys/devices/system/cpu/cpu{CPU}/cpufreq/scaling_governor') as f:
            gov = f.read().strip()
    except OSError:
        gov = 'unknown'
    return {'commit': commit, 'load': os.getloadavg()[0], 'governor': gov, 'cpu': CPU}


def paired_summary(parent, child, n_boot=4000, seed=0):
    """Per-round PAIRED ratios child/parent: `(median, lo, hi, wins)` --
    the median ratio, a bootstrap 95 % interval of that median over the
    rounds, and how many rounds the child was faster.  Pairing per round
    cancels the drift a pooled comparison of two sequential runs
    attributes to the code (pyperf's own significance test is unpaired)."""
    import random
    import statistics
    r = [c / p for p, c in zip(parent, child)]
    if not r:
        return None
    med = statistics.median(r)
    rng = random.Random(seed)
    boots = sorted(statistics.median(rng.choices(r, k=len(r))) for _ in range(n_boot))
    lo, hi = boots[int(0.025 * n_boot)], boots[int(0.975 * n_boot) - 1]
    return med, lo, hi, sum(1 for v in r if v < 1.0)

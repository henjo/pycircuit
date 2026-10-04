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


def equal_paths(trees):
    """`trees` (`{label: path}`) reached through symlinks of ONE length
    (`<tmp>/p`, `<tmp>/c`), removed at exit.  A module's path is a string
    the process allocates, so two trees at paths of different lengths lay
    their heaps out differently, and numpy's identity hashing and the
    object-keyed dicts probe differently: the CPU counter read -0.56 to
    +0.10 % between two trees at ONE commit, -0.08 to +0.23 % through
    equal-length links (2026-10-04)."""
    import atexit
    import shutil
    import tempfile
    d = tempfile.mkdtemp(prefix='pyc-trees-')
    atexit.register(shutil.rmtree, d, True)
    out = {}
    for lab, path in trees.items():
        link = os.path.join(d, lab[0])
        os.symlink(os.path.abspath(path), link)
        out[lab] = link
    return out


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


## -- instruction counts (2026-10-04, testing for development, stage 6) ---------------
##
## Wall time on this shared box resolves about 1 % with eight paired rounds;
## the instructions a call retires are nearly a property of the code alone
## (with address randomisation off and a fixed hash seed, three cachegrind
## runs gave the same count to the instruction).  Two routes:
##
## * `perf`: the CPU's own counter, read around the call in-process
##   (`perf_event_open`, user space only, the P-core's PMU on this hybrid
##   CPU).  Needs `kernel.perf_event_paranoid` <= 2 -- one root setting:
##       sudo sysctl kernel.perf_event_paranoid=2
##       echo kernel.perf_event_paranoid=2 | sudo tee /etc/sysctl.d/60-perf.conf
## * `valgrind`: cachegrind, its counting switched on only around the call
##   (`scripts/get_valgrind.sh` unpacks it without root).  Every process runs
##   under valgrind's translator: minutes per process, so one process per case.
##
## `count_route()` names the route this box offers (perf first) or None.

PERF_TYPE_HARDWARE = 0
PERF_COUNT_HW_INSTRUCTIONS = 1
PERF_FORMAT_TOTAL_TIME_ENABLED = 1
PERF_FORMAT_TOTAL_TIME_RUNNING = 2
PERF_EVENT_IOC_ENABLE = 0x2400
PERF_EVENT_IOC_DISABLE = 0x2401
PERF_EVENT_IOC_RESET = 0x2403
_NR_PERF_EVENT_OPEN = 298                      # x86_64
VALGRIND_DIR = os.environ.get('PYCIRCUIT_VALGRIND') or os.path.join(
    os.path.expanduser('~'), '.local', 'opt', 'valgrind')


def _paranoid():
    try:
        with open('/proc/sys/kernel/perf_event_paranoid') as f:
            return int(f.read().strip())
    except (OSError, ValueError):
        return None


def _core_pmu_type():
    """The P-core PMU's type on a hybrid CPU (`cpu_core`), or None."""
    try:
        with open('/sys/bus/event_source/devices/cpu_core/type') as f:
            return int(f.read().strip())
    except (OSError, ValueError):
        return None


def _perf_attr(config):
    """`struct perf_event_attr`, version 0 (64 bytes): type, size, config,
    sample_period, sample_type, read_format (time enabled and running), the
    flag bits (disabled, exclude_kernel, exclude_hv), wakeup_events, bp_type,
    config1."""
    import struct
    flags = 1 | (1 << 5) | (1 << 6)
    return struct.pack('IIQQQQQIIQ', PERF_TYPE_HARDWARE, 64, config, 0, 0,
                       PERF_FORMAT_TOTAL_TIME_ENABLED | PERF_FORMAT_TOTAL_TIME_RUNNING,
                       flags, 0, 0, 0)


class InstrCounter:
    """User-space instructions retired by this process between `start()`
    and `stop()`, from the CPU's counter.  `count` is None when the counter
    did not run the whole time (multiplexed, or the process left the P-core
    whose PMU it counts) -- such a reading is rejected, not scaled."""

    def __init__(self):
        import ctypes
        import struct
        self._ct, self._st = ctypes, struct
        libc = ctypes.CDLL(None, use_errno=True)
        self._libc = libc
        config = PERF_COUNT_HW_INSTRUCTIONS
        pmu = _core_pmu_type()
        if pmu is not None:
            config |= pmu << 32               # (the extended type: one PMU of a hybrid CPU)
        attr = _perf_attr(config)
        assert len(attr) == 64
        buf = ctypes.create_string_buffer(attr, 64)
        libc.syscall.restype = ctypes.c_long
        fd = libc.syscall(ctypes.c_long(_NR_PERF_EVENT_OPEN), buf, ctypes.c_int(0),
                          ctypes.c_int(-1), ctypes.c_int(-1), ctypes.c_ulong(8))
        if fd < 0:
            err = ctypes.get_errno()
            raise OSError(err, f'perf_event_open: {os.strerror(err)} '
                               f'(perf_event_paranoid={_paranoid()})')
        self.fd = fd
        self.count = None

    def start(self):
        self._libc.ioctl(self.fd, PERF_EVENT_IOC_RESET, 0)
        self._libc.ioctl(self.fd, PERF_EVENT_IOC_ENABLE, 0)

    def stop(self):
        self._libc.ioctl(self.fd, PERF_EVENT_IOC_DISABLE, 0)
        value, enabled, running = self._st.unpack('QQQ', os.read(self.fd, 24))
        self.count = value if running == enabled and enabled > 0 else None
        return self.count

    def close(self):
        os.close(self.fd)


def perf_available():
    try:
        InstrCounter().close()
    except OSError:
        return False
    return True


def valgrind_available():
    return (os.path.exists(os.path.join(VALGRIND_DIR, 'usr', 'bin', 'valgrind.bin'))
            and os.path.exists(os.path.join(VALGRIND_DIR, 'pycircuit-count.so')))


def count_route():
    """'perf', 'valgrind' or None (and why)."""
    if perf_available():
        return 'perf'
    if valgrind_available():
        return 'valgrind'
    return None


def count_unavailable_message():
    return ('no instruction counter here: the kernel keeps its counters closed '
            f'(perf_event_paranoid={_paranoid()}; `sudo sysctl kernel.perf_event_paranoid=2` '
            'opens them for your own processes) and no valgrind in '
            f'{VALGRIND_DIR} (`scripts/get_valgrind.sh` unpacks it, no root)')


def valgrind_prefix(out_file):
    """The command prefix running a child under cachegrind, counting only
    between the helper's start and stop (`--instr-at-start=no`)."""
    return [os.path.join(VALGRIND_DIR, 'usr', 'bin', 'valgrind.bin'), '--tool=cachegrind',
            '--cache-sim=no', '--instr-at-start=no', f'--cachegrind-out-file={out_file}']


def valgrind_env(env):
    env = dict(env)
    env['VALGRIND_LIB'] = os.path.join(VALGRIND_DIR, 'usr', 'libexec', 'valgrind')
    return env


def valgrind_count(out_file):
    """The instructions cachegrind counted (its `summary:` line)."""
    with open(out_file) as f:
        for line in f:
            if line.startswith('summary:'):
                return int(line.split()[1])
    raise ValueError(f'no summary in {out_file}')


def split_counts(c1, c2, steps):
    """`(marginal step, per-solve cost)` from the instructions of one
    circuit's `steps`-step solve `c1` and its `2 * steps`-step solve `c2`:
    the second solve's extra steps cost `c2 - c1`, the rest of `c1` is what
    a solve pays once (its setup, the operating point, the result).  None
    where a count is missing (speed round 8, stage 1)."""
    if c1 is None or c2 is None:
        return None
    step = (c2 - c1) / steps
    return step, c1 - steps * step


class CountRegion:
    """In a child: the region whose instructions are counted, by the route
    the parent chose (`PYCIRCUIT_BENCH_COUNT`): `perf` reads the counter
    (`.count`), `valgrind` switches cachegrind on and off (the parent reads
    the file).  Outside a counting run it does nothing."""

    def __init__(self):
        self.route = os.environ.get('PYCIRCUIT_BENCH_COUNT') or None
        self.count = None
        self._perf = self._helper = None
        if self.route == 'perf':
            self._perf = InstrCounter()
        elif self.route == 'valgrind':
            import ctypes
            self._helper = ctypes.CDLL(os.path.join(VALGRIND_DIR, 'pycircuit-count.so'))

    def __enter__(self):
        if self._perf is not None:
            self._perf.start()
        elif self._helper is not None:
            self._helper.pyc_count_start()
        return self

    def __exit__(self, *exc):
        if self._perf is not None:
            self.count = self._perf.stop()
        elif self._helper is not None:
            self._helper.pyc_count_stop()
        return False

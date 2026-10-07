"""The public benchmark circuits, in the local cache -- never in the repository.

The CircuitSim90/MCNC circuits come from Sandia's Xyce regression suite, which
declares no license; the IBM power grids' page states no terms.  So neither is
vendored: `benchmarks/fetch_spice_suite.py` downloads them into the cache,
pinned and hashed by the manifest beside this module
(`spice_suite_manifest.json`: every file's source, path, size and sha256, each
source's pinned version and URL).  Everything that reads the data goes through
`spice_data`, which answers None where a file is absent, so a test or a
benchmark skips instead of reaching for the network.

The cache: ``$PYCIRCUIT_BENCH_DATA``, else ``$XDG_CACHE_HOME/pycircuit/
benchmarks``, else ``~/.cache/pycircuit/benchmarks``; a file lives at
``<cache>/<source>/<version>/<path>``.
"""
import functools
import json
import os

MANIFEST = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'spice_suite_manifest.json')


def cache_root():
    """Where the fetched files live (see the module note)."""
    env = os.environ.get('PYCIRCUIT_BENCH_DATA')
    if env:
        return env
    xdg = os.environ.get('XDG_CACHE_HOME') or os.path.join(os.path.expanduser('~'), '.cache')
    return os.path.join(xdg, 'pycircuit', 'benchmarks')


@functools.cache
def manifest():
    """The manifest, read once: ``{'sources': {...}, 'files': [...]}``."""
    with open(MANIFEST, encoding='utf-8') as fh:
        return json.load(fh)


def entry(path, source='xyce_regression'):
    """The manifest's record of `path` in `source`, or None."""
    for f in manifest()['files']:
        if f['source'] == source and f['path'] == path:
            return f
    return None


def local_path(path, source='xyce_regression'):
    """Where `path` of `source` lives in the cache (whether or not it is there)."""
    version = manifest()['sources'][source]['version']
    return os.path.join(cache_root(), source, version, *path.split('/'))


def spice_data(path, source='xyce_regression'):
    """The cached file `path` of `source`, or None where it is absent or not
    the manifest's size (a partial download) -- the one gate every
    data-dependent test and benchmark goes through.  A path the manifest does
    not list is a ValueError: a typo should not read as missing data."""
    rec = entry(path, source)
    if rec is None:
        raise ValueError(f'{source}:{path} is not in the benchmark manifest')
    p = local_path(path, source)
    try:
        if os.path.getsize(p) != rec['size']:
            return None
    except OSError:
        return None
    return p

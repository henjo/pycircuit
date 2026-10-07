"""`pycircuit._testing.benchdata` and `benchmarks/fetch_spice_suite.py`: the
public benchmark circuits live in a local cache, pinned by a committed
manifest, never in the repository.  The manifest is well formed (every
source pinned to a version, every file a relative path with its size and
sha256); `spice_data` answers a cached file of the manifest's size, None for
an absent or partial one, and a ValueError for a path the manifest does not
list; the fetcher stores only the manifest's bytes.  No test here touches
the network (the fetcher's download is replaced)."""
import hashlib
import importlib.util
import os
import re

import pytest

from pycircuit._testing import benchdata

ROOT = os.path.dirname(os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__)))))
FETCH = os.path.join(ROOT, 'benchmarks', 'fetch_spice_suite.py')


def test_the_manifest_is_well_formed():
    man = benchdata.manifest()
    assert man['format'] == 1 and set(man['sources']) == {'xyce_regression', 'ibm_pg'}
    for name, src in man['sources'].items():
        assert {'what', 'home', 'version', 'url', 'terms'} <= set(src), name
        assert '{path}' in src['url'], name
        assert src['url'].startswith('https://') and 'never vendored' in src['terms'], name
    ## (the Xyce files are fetched AT a commit, not a branch: the bytes its
    ## version names cannot move; the IBM page has no versions -- its
    ## version a date label, the sha256 the only pin: moved bytes are refused)
    xyce = man['sources']['xyce_regression']
    assert re.fullmatch('[0-9a-f]{40}', xyce['version']) and '{version}' in xyce['url']
    seen = set()
    for f in man['files']:
        key = (f['source'], f['path'])
        assert key not in seen, key
        seen.add(key)
        assert f['source'] in man['sources'], key
        parts = f['path'].split('/')
        assert f['path'] == f['path'].strip() and '\\' not in f['path'], key
        assert parts[0] and '..' not in parts and '.' not in parts and '' not in parts, key
        assert isinstance(f['size'], int) and f['size'] > 0, key
        assert re.fullmatch('[0-9a-f]{64}', f['sha256']), key
    counts = {}
    for f in man['files']:
        k = f['source'] if f['source'] != 'xyce_regression' else f['path'].split('/')[0]
        counts[k] = counts.get(k, 0) + 1
    assert counts == {'Netlists': 52, 'OutputData': 30, 'ibm_pg': 2}, counts


def test_no_benchmark_file_is_in_the_repository():
    """No file of the package, the benchmarks or the documents carries a
    manifest file's bytes (sizes first, a hash where one matches)."""
    by_size = {}
    for f in benchdata.manifest()['files']:
        by_size.setdefault(f['size'], set()).add(f['sha256'])
    found = []
    for top in ('pycircuit', 'benchmarks', 'doc', 'scripts'):
        for dirpath, dirnames, filenames in os.walk(os.path.join(ROOT, top)):
            dirnames[:] = [d for d in dirnames if d != '__pycache__']
            for fn in filenames:
                p = os.path.join(dirpath, fn)
                try:
                    size = os.path.getsize(p)
                except OSError:
                    continue
                if size in by_size:
                    with open(p, 'rb') as fh:
                        if hashlib.sha256(fh.read()).hexdigest() in by_size[size]:
                            found.append(os.path.relpath(p, ROOT))
    assert not found, found


def test_the_cache_root_follows_the_environment(monkeypatch, tmp_path):
    monkeypatch.setenv('PYCIRCUIT_BENCH_DATA', str(tmp_path / 'own'))
    monkeypatch.setenv('XDG_CACHE_HOME', str(tmp_path / 'xdg'))
    assert benchdata.cache_root() == str(tmp_path / 'own')
    monkeypatch.delenv('PYCIRCUIT_BENCH_DATA')
    assert benchdata.cache_root() == str(tmp_path / 'xdg' / 'pycircuit' / 'benchmarks')
    monkeypatch.delenv('XDG_CACHE_HOME')
    monkeypatch.setenv('HOME', str(tmp_path / 'home'))
    assert benchdata.cache_root() == str(tmp_path / 'home' / '.cache' / 'pycircuit' / 'benchmarks')


def _first(source):
    return next(f for f in benchdata.manifest()['files'] if f['source'] == source)


@pytest.mark.parametrize('source', ['xyce_regression', 'ibm_pg'])
def test_spice_data_answers_the_cached_file_or_none(monkeypatch, tmp_path, source):
    """Absent: None.  Shorter or longer than the manifest's size (a partial
    download): None.  The manifest's size: the path, under
    ``<cache>/<source>/<version>/<path>``."""
    monkeypatch.setenv('PYCIRCUIT_BENCH_DATA', str(tmp_path))
    f = _first(source)
    version = benchdata.manifest()['sources'][source]['version']
    want = os.path.join(str(tmp_path), source, version, *f['path'].split('/'))
    assert benchdata.local_path(f['path'], source) == want
    assert benchdata.spice_data(f['path'], source) is None
    os.makedirs(os.path.dirname(want))
    for size in (f['size'] - 1, f['size'] + 1, 0):
        with open(want, 'wb') as fh:
            fh.write(b'x' * size)
        assert benchdata.spice_data(f['path'], source) is None, size
    with open(want, 'wb') as fh:
        fh.write(b'x' * f['size'])
    assert benchdata.spice_data(f['path'], source) == want


def test_a_path_the_manifest_does_not_list_is_an_error():
    """A typo must not read as missing data (a skip)."""
    with pytest.raises(ValueError, match='not in the benchmark manifest'):
        benchdata.spice_data('Netlists/CircuitSim90/MOS2/no_such.cir')
    f = _first('ibm_pg')
    with pytest.raises(ValueError, match='not in the benchmark manifest'):
        benchdata.spice_data(f['path'])           # the right path, the wrong source


def _fetcher():
    spec = importlib.util.spec_from_file_location('_fetch_spice_suite', FETCH)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


def test_the_fetcher_stores_only_the_manifests_bytes(monkeypatch, tmp_path, capsys):
    """Every download is checked against the manifest's sha256 and size: the
    wrong bytes are refused (nothing stored, the exit status 1); the right
    ones land at `local_path`, atomically (no temporary file left); a cached
    file of the right size is not fetched again; `verify` re-fetches a file
    whose hash differs.  The download is replaced: no network."""
    monkeypatch.setenv('PYCIRCUIT_BENCH_DATA', str(tmp_path))
    fs = _fetcher()
    man = benchdata.manifest()
    picked = [_first('xyce_regression'), _first('ibm_pg')]
    ## (the manifest's hashes stand for bytes these tests do not have: each
    ## picked file's bytes are made up and the manifest pinned to them)
    data = {f['path']: os.urandom(f['size']) for f in picked}
    picked = [{**f, 'sha256': hashlib.sha256(data[f['path']]).hexdigest()} for f in picked]
    monkeypatch.setattr(benchdata, 'manifest', lambda: {**man, 'files': picked})
    asked = []

    def wrong(url, timeout=120):
        ## (the manifest's size, other bytes: only the hash can refuse them)
        asked.append(url)
        return bytes(next(f['size'] for f in picked if url.endswith('/' + f['path'])))
    monkeypatch.setattr(fs, '_get', wrong)
    assert fs.fetch() == 1 and len(asked) == 2
    assert 'REFUSED' in capsys.readouterr().out
    assert not any(os.path.exists(benchdata.local_path(f['path'], f['source'])) for f in picked)

    def right(url, timeout=120):
        asked.append(url)
        return next(data[f['path']] for f in picked if url.endswith('/' + f['path']))
    asked.clear()
    monkeypatch.setattr(fs, '_get', right)
    assert fs.fetch() == 0 and len(asked) == 2
    for f in picked:
        src = man['sources'][f['source']]
        assert src['url'].format(version=src['version'], path=f['path']) in asked
        with open(benchdata.local_path(f['path'], f['source']), 'rb') as fh:
            assert fh.read() == data[f['path']]
    left = [fn for _d, _s, fns in os.walk(tmp_path) for fn in fns if fn.startswith('.part-')]
    assert not left, left
    asked.clear()
    assert fs.fetch() == 0 and not asked                   # all cached: nothing fetched
    p = benchdata.local_path(picked[0]['path'], picked[0]['source'])
    with open(p, 'r+b') as fh:                              # same size, other bytes
        fh.write(b'\x00' if data[picked[0]['path']][:1] != b'\x00' else b'\x01')
    assert fs.fetch() == 0 and not asked                   # a size check cannot see it
    assert fs.fetch(verify=True) == 0 and len(asked) == 1  # a hash can
    with open(p, 'rb') as fh:
        assert fh.read() == data[picked[0]['path']]

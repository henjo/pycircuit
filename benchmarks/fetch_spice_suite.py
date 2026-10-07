"""Fetch the public SPICE benchmark circuits into the local cache (never into
the repository): the CircuitSim90 and MCNC circuits with Xyce's gold outputs
from Sandia's Xyce regression suite, two of its harmonic-balance circuits, the
4049 oscillator, and the smallest IBM transient power grid with its reference
output.  Every file is pinned by `pycircuit/_testing/spice_suite_manifest.json`
(its source's version, its size and sha256); `pycircuit._testing.benchdata`
reads them.

    python benchmarks/fetch_spice_suite.py            # fetch what is missing
    python benchmarks/fetch_spice_suite.py --verify   # re-hash every cached file
    python benchmarks/fetch_spice_suite.py --build-manifest   # maintainers: list,
                                                      # download and hash anew

Why not vendored: the Xyce regression suite declares no license and the IBM
page states no terms; a URL and a hash are ours to keep, the files are not.
"""
import hashlib
import json
import os
import re
import sys
import tempfile
import urllib.request

if __name__ == '__main__':
    ## (run as a script: the repository root importable; imported by the
    ## tests it is already, and an import must not change `sys.path`)
    sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from pycircuit._testing import benchdata

XYCE_VERSION = '7bb7e98f0ed3a81a7d1cf1d10b68592107ed40b2'
SOURCES = {
    'xyce_regression': {
        'what': "Sandia's Xyce regression suite: CircuitSim90 / MCNC circuits and Xyce gold outputs",
        'home': 'https://github.com/Xyce/Xyce_Regression',
        'version': XYCE_VERSION,
        'url': 'https://raw.githubusercontent.com/Xyce/Xyce_Regression/{version}/{path}',
        'terms': 'no license declared: fetched into the cache, never vendored',
    },
    'ibm_pg': {
        'what': 'IBM power grid benchmarks (S. Nassif), the transient set',
        'home': 'https://web.ece.ucsb.edu/~lip/PGBenchmarks/ibmpgbench.html',
        'version': '2020-03',
        'url': 'https://web.ece.ucsb.edu/~lip/PGBenchmarks/ibmpg/{path}',
        'terms': 'no terms stated: fetched into the cache, never vendored',
    },
}
#: which files of the Xyce suite (paths at `XYCE_VERSION`)
XYCE_SELECT = [
    r'Netlists/CircuitSim90/(BJT|MOS2|MOS3|MOS2_LARGE)/[^/]+\.cir',
    r'Netlists/CircuitSim90/MOS2_LARGE/(README|voter\.cir\.options)',
    r'Netlists/MCNC_BJT_(LATCH|OPAMPAL|RCA|SCHMITECL)/[^/]+\.cir(_NORUN)?',
    r'Netlists/4049OSC/4049osc\.cir',
    r'Netlists/HB/(common_emitter_hb|gilbert_cell_hb)\.cir',
    r'OutputData/CircuitSim90/(BJT|MOS2|MOS3|MOS2_LARGE)/[^/]+\.prn',
    r'OutputData/MCNC_BJT_(LATCH|OPAMPAL|RCA|SCHMITECL)/[^/]+\.prn',
    r'OutputData/4049OSC/[^/]+\.prn',
    r'OutputData/HB/(common_emitter_hb|gilbert_cell_hb)\.cir\.[^/]+\.prn',
]
IBM_FILES = ['ibmpg1t.spice.bz2', 'ibmpg1t.output.bz2']


def _get(url, timeout=120):
    req = urllib.request.Request(url, headers={'User-Agent': 'pycircuit-benchmark-fetch'})
    with urllib.request.urlopen(req, timeout=timeout) as r:
        return r.read()


def _store(path, data):
    """`data` at `path`, atomically (a temporary file in the same directory
    moved into place): a reader never sees a partial file."""
    os.makedirs(os.path.dirname(path), exist_ok=True)
    fd, tmp = tempfile.mkstemp(dir=os.path.dirname(path), prefix='.part-')
    try:
        with os.fdopen(fd, 'wb') as fh:
            fh.write(data)
        os.replace(tmp, path)
    except BaseException:
        if os.path.exists(tmp):
            os.unlink(tmp)
        raise


def _xyce_paths():
    tree = json.loads(_get('https://api.github.com/repos/Xyce/Xyce_Regression/git/trees/'
                           f'{XYCE_VERSION}?recursive=1'))
    if tree.get('truncated'):
        raise RuntimeError('the GitHub tree listing came back truncated')
    pats = [re.compile(p + r'\Z') for p in XYCE_SELECT]
    return sorted(e['path'] for e in tree['tree']
                  if e['type'] == 'blob' and any(p.match(e['path']) for p in pats))


def build_manifest():
    """List, download and hash every selected file; write the manifest (the
    files land in the cache as they are hashed)."""
    files, got = [], []
    wanted = [('xyce_regression', p) for p in _xyce_paths()] + [('ibm_pg', p) for p in IBM_FILES]
    for k, (source, path) in enumerate(wanted, 1):
        data = _get(SOURCES[source]['url'].format(version=SOURCES[source]['version'], path=path))
        files.append({'source': source, 'path': path, 'size': len(data),
                      'sha256': hashlib.sha256(data).hexdigest()})
        got.append(data)
        print(f'  [{k}/{len(wanted)}] {source}:{path} {len(data)} bytes', flush=True)
    out = {'format': 1, 'sources': SOURCES, 'files': files}
    with open(benchdata.MANIFEST, 'w', encoding='utf-8') as fh:
        json.dump(out, fh, indent=1, sort_keys=True)
        fh.write('\n')
    benchdata.manifest.cache_clear()
    ## (stored once the manifest names the version they live under)
    for f, data in zip(files, got, strict=True):
        _store(benchdata.local_path(f['path'], f['source']), data)
    print(f'manifest: {len(files)} files, {sum(f["size"] for f in files)} bytes')


def fetch(verify=False):
    """Fetch every manifest file the cache lacks (or, `verify`, whose hash
    differs); a download whose hash is not the manifest's is refused."""
    man = benchdata.manifest()
    got = bad = 0
    for f in man['files']:
        src = man['sources'][f['source']]
        path = benchdata.local_path(f['path'], f['source'])
        if os.path.exists(path):
            if not verify:
                if os.path.getsize(path) == f['size']:
                    continue
            else:
                with open(path, 'rb') as fh:
                    if hashlib.sha256(fh.read()).hexdigest() == f['sha256']:
                        continue
                print(f'  hash differs, fetching again: {f["source"]}:{f["path"]}')
        data = _get(src['url'].format(version=src['version'], path=f['path']))
        if hashlib.sha256(data).hexdigest() != f['sha256'] or len(data) != f['size']:
            print(f'  REFUSED (not the manifest\'s bytes): {f["source"]}:{f["path"]}')
            bad += 1
            continue
        _store(path, data)
        got += 1
    print(f'{got} fetched, {bad} refused, {len(man["files"])} in the manifest; '
          f'cache {benchdata.cache_root()}')
    return 1 if bad else 0


if __name__ == '__main__':
    if '--build-manifest' in sys.argv:
        build_manifest()
        sys.exit(0)
    sys.exit(fetch(verify='--verify' in sys.argv))

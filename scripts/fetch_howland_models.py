#!/usr/bin/env python
"""Download the vendor SPICE models examples 12 and 13 need.

The models' licences do not allow committing them to this repository, so
they are fetched from the vendors instead, into
``doc/src/circuit/examples/howland/`` (or ``$PYCIRCUIT_HOWLAND_MODELS``):

``PBSS4230PAN.txt``
    Nexperia's model, a direct download.
``PC817.sub``
    Not published on its own: it is in the ``lib.zip`` inside LTspice's
    Windows installer.  The installer (about 145 MB) is an MSI -- an OLE
    compound file holding an MSZIP cabinet -- and is unpacked here with
    `olefile` and `zlib`, so no Windows and no system tools are needed.
``ad8606.cir``
    Analog Devices' macromodel.  www.analog.com does not answer scripted
    requests from every network (seen 2026-10-06: every request hangs, the
    home page included), so this one may fail.  That is a WARNING, not an
    error: without it the documentation shows the results recorded when the
    model was present, and says so on the page.

Every file is checked for the subcircuit it must define -- a vendor site
that answers with an HTML error page must not pass as a model -- and
written as UTF-8 with LF line ends, which is how `spicecard` reads it.

Files already present are kept unless ``--force``.  The extracted
``PC817.sub`` is also kept in ``--cache`` (default ``~/.cache/pycircuit``)
so the installer is downloaded once, not on every build.

Exit status is 0 when everything that could be fetched was, even if the
AD8606 was not; ``--strict`` makes any missing model an error.
"""
import argparse
import io
import os
import struct
import sys
import urllib.request
import zipfile
import zlib

HERE = os.path.dirname(os.path.abspath(__file__))
DEST = os.environ.get(
    'PYCIRCUIT_HOWLAND_MODELS',
    os.path.join(HERE, os.pardir, 'doc', 'src', 'circuit', 'examples',
                 'howland'))

NEXPERIA = 'https://assets.nexperia.com/documents/spice-model/PBSS4230PAN.txt'
ADI = 'https://www.analog.com/media/en/simulation-models/spice-models/ad8606.cir'
LTSPICE_MSI = 'https://ltspice.analog.com/software/LTspice64.msi'

#: What each file must contain to be the model and not an error page.
MUST_DEFINE = {'PBSS4230PAN.txt': '.subckt pbss4230pan',
               'ad8606.cir': '.subckt ad8606',
               'PC817.sub': '.subckt pc817'}

UA = ('Mozilla/5.0 (X11; Linux x86_64) AppleWebKit/537.36 '
      '(KHTML, like Gecko) Chrome/130.0 Safari/537.36')


class FetchError(Exception):
    pass


def _get(url, timeout):
    req = urllib.request.Request(url, headers={'User-Agent': UA,
                                               'Accept': '*/*'})
    try:
        with urllib.request.urlopen(req, timeout=timeout) as r:
            return r.read()
    except Exception as exc:
        raise FetchError('%s: %s' % (url, exc)) from exc


def _normalise(name, raw):
    """UTF-8, LF line ends, and the subcircuit it must define."""
    try:
        text = raw.decode('utf-8')
    except UnicodeDecodeError:
        text = raw.decode('latin-1')         # LTspice's library files
    text = text.replace('\r\n', '\n')
    if MUST_DEFINE[name] not in text.lower():
        raise FetchError('%s does not define %r -- not a model file '
                         '(an error page?)' % (name, MUST_DEFINE[name]))
    return text.encode('utf-8')


## ------------------------------------------------------------ LTspice MSI
def _cab_file(cab, wanted):
    """One file out of an MSZIP cabinet, by name."""
    (_, _, _, _, coff_files, _, _, _, n_folders, n_files,
     flags) = struct.unpack('<4sIIIIIBBHHH', cab[:32])
    off = 36 + (4 if flags & 4 else 0)
    folders = []
    for _ in range(n_folders):
        folders.append(struct.unpack('<IHH', cab[off:off + 8]))
        off += 8
    p, files = coff_files, []
    for _ in range(n_files):
        size, uoff, ifolder, _, _, _ = struct.unpack('<IIHHHH',
                                                     cab[p:p + 16])
        end = cab.index(b'\0', p + 16)
        files.append((cab[p + 16:end].decode('latin-1'), size, uoff,
                      ifolder))
        p = end + 1
    hits = [f for f in files if f[0] == wanted]
    if not hits:
        raise FetchError('no %r in the LTspice installer' % wanted)
    _, size, uoff, ifolder = hits[0]
    start, n_blocks, compress = folders[ifolder]
    if compress & 0xF != 1:
        raise FetchError('cabinet is not MSZIP (type %d)' % compress)
    ## MSZIP: each block is 'CK' + raw deflate, with the previous block's
    ## output as its dictionary.
    out, q, prev = bytearray(), start, b''
    for _ in range(n_blocks):
        _, cb, _ = struct.unpack('<IHH', cab[q:q + 8])
        block = cab[q + 10:q + 8 + cb]
        q += 8 + cb
        d = (zlib.decompressobj(-15, zdict=prev) if prev
             else zlib.decompressobj(-15))
        prev = d.decompress(block) + d.flush()
        out += prev
        if len(out) >= uoff + size:
            break
    return bytes(out[uoff:uoff + size])


def pc817_from_ltspice(timeout):
    try:
        import olefile
    except ImportError as exc:
        raise FetchError('reading the LTspice installer needs olefile '
                         '(pip install olefile)') from exc
    print('downloading the LTspice installer (about 145 MB) ...',
          flush=True)
    msi = olefile.OleFileIO(io.BytesIO(_get(LTSPICE_MSI, timeout)))
    ## The cabinet is the MSI's largest stream.
    cab = msi.openstream(max(msi.listdir(), key=msi.get_size)).read()
    lib = zipfile.ZipFile(io.BytesIO(_cab_file(cab, 'lib.zip')))
    return lib.read('lib/sub/PC817.sub')


## ------------------------------------------------------------------- main
def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__.split('\n')[0])
    ap.add_argument('--dest', default=DEST)
    ap.add_argument('--cache', default=os.path.join(
        os.path.expanduser('~'), '.cache', 'pycircuit'))
    ap.add_argument('--force', action='store_true',
                    help='download again even if the file is present')
    ap.add_argument('--strict', action='store_true',
                    help='fail if any model could not be fetched')
    ap.add_argument('--timeout', type=float, default=60.0)
    a = ap.parse_args(argv)
    os.makedirs(a.dest, exist_ok=True)
    os.makedirs(a.cache, exist_ok=True)

    def nexperia():
        return _get(NEXPERIA, a.timeout)

    def adi():
        return _get(ADI, a.timeout)

    def pc817():
        cached = os.path.join(a.cache, 'PC817.sub')
        if os.path.isfile(cached) and not a.force:
            return open(cached, 'rb').read()
        raw = pc817_from_ltspice(max(a.timeout, 600.0))
        open(cached, 'wb').write(raw)
        return raw

    missing = []
    for name, fetch in (('PBSS4230PAN.txt', nexperia),
                        ('PC817.sub', pc817), ('ad8606.cir', adi)):
        path = os.path.join(a.dest, name)
        if os.path.isfile(path) and not a.force:
            print('%-16s present' % name)
            continue
        try:
            data = _normalise(name, fetch())
        except FetchError as exc:
            missing.append(name)
            msg = 'could not fetch %s: %s' % (name, exc)
            ## A GitHub Actions annotation, and plain text elsewhere.
            print(('::warning::' if os.environ.get('GITHUB_ACTIONS')
                   else 'WARNING: ') + msg, file=sys.stderr)
            continue
        with open(path, 'wb') as f:
            f.write(data)
        print('%-16s fetched' % name)

    if missing:
        print('missing: %s -- examples 12 and 13 will show their recorded '
              'results' % ', '.join(missing), file=sys.stderr)
    return 1 if (missing and a.strict) else 0


if __name__ == '__main__':
    sys.exit(main())

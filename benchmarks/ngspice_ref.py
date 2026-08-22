"""Run an ngspice deck and bring its data back -- the golden-reference harness.

pycircuit's compact-model work is checked against real models: the IHP
Open PDK ships PSP103, r3_cmc, mosvar and the two MoM capacitors as
Verilog-A plus compiled OSDI binaries, so each can be swept in ngspice
and the curves kept as a fixed artifact in the repo.

The tests never call this.  A reference is generated once, committed, and
read from disk, so a test failure means pycircuit moved rather than
someone's PDK checkout.

Two traps this exists to encapsulate, both of which produced a reference
that loaded fine and was wrong:

* ngspice's ``print`` SPLITS into separate tables once the columns exceed
  the terminal width, so a parser reading the first table silently drops
  the trailing vectors.  Everything here goes through ``wrdata``, which
  writes one row per point regardless of how many vectors there are.
* ``wrdata`` emits a ``(scale, value)`` PAIR per vector, and for a
  complex vector (any ``ac`` result) a ``(scale, real, imag)`` TRIPLE.
  Column counting has to know which.
"""

import os
import subprocess
import tempfile


DEFAULT_PDK = os.path.expanduser(
    '~/source/IHP-Open-PDK/ihp-sg13g2/libs.tech/ngspice')


def version():
    """The ngspice banner line, recorded in every reference file."""
    out = subprocess.run(['ngspice', '--version'], capture_output=True,
                         text=True).stdout.splitlines()
    return next((l.strip() for l in out if 'ngspice-' in l), 'unknown')


def run(deck, nvectors, complex_data=False, timeout=300, label=''):
    """Run `deck` (a format string wanting `{out}`) and return its columns.

    `deck` must contain a ``wrdata {out} <vectors...>`` line listing
    exactly `nvectors` vectors.  Returns a list of `nvectors` lists, plus
    the sweep column, as ``(sweep, [v0, v1, ...])``.  With
    `complex_data`, each vector comes back as ``(real, imag)`` pairs.
    """
    tmpdir = tempfile.mkdtemp()
    datafile = os.path.join(tmpdir, 'out.dat')
    deckfile = os.path.join(tmpdir, 'deck.sp')
    try:
        with open(deckfile, 'w') as fh:
            fh.write(deck.format(out=datafile))
        proc = subprocess.run(['ngspice', '-b', deckfile],
                              capture_output=True, text=True,
                              timeout=timeout)
        if proc.returncode != 0:
            raise RuntimeError('ngspice failed for %s:\n%s'
                               % (label, proc.stderr[-2000:]))
        if not os.path.exists(datafile):
            raise RuntimeError('ngspice wrote no data for %s:\n%s'
                               % (label, proc.stdout[-2000:]))

        ## `wrdata` writes (scale, value) per real vector and
        ## (scale, real, imag) per complex one.
        per = 3 if complex_data else 2
        want = nvectors * per
        sweep, cols = [], [[] for _ in range(nvectors)]
        with open(datafile) as fh:
            for line in fh:
                parts = line.split()
                if len(parts) != want:
                    continue
                vals = [float(v) for v in parts]
                sweep.append(vals[0])
                for k in range(nvectors):
                    chunk = vals[k * per:(k + 1) * per]
                    cols[k].append(tuple(chunk[1:]) if complex_data
                                   else chunk[1])
        if not sweep:
            raise RuntimeError(
                'no data rows parsed for %s: expected %d columns per row'
                % (label, want))
        return sweep, cols
    finally:
        for f in (datafile, deckfile):
            if os.path.exists(f):
                os.unlink(f)
        if os.path.isdir(tmpdir):
            os.rmdir(tmpdir)

"""Read circuit simulators' outputs, and compare waveforms as Xyce's
regression suite compares them (the SPICE benchmark plan's stage 4).

Readers:

- `read_prn`: Xyce's ``.prn`` -- a header line (``Index TIME col ...``), a
  row per output point, ``End of Xyce(TM) Simulation``; the harmonic
  balance's ``.HB.TD.prn`` and ``.HB.FD.prn`` (``TIME`` / ``FREQ``,
  ``Re(...)`` / ``Im(...)``) have the same layout.  Returns the column
  names and a float array.
- `read_ngspice_print`: ngspice's batch ``.print`` output -- tables of
  ``Index time v(a) ...``, the columns split over tables by the page
  width, repeated per page.  Returns the column names and a float array,
  one row per index.
- `read_ibm_output`: the IBM power grid benchmarks' reference (``Node:
  name`` blocks of ``time volts`` pairs): {node: (t, v)}.

`xyce_verify` is the transient comparison of ``xyce_verify.pl`` (Sandia's
Xyce regression suite, ``TestScripts/xyce_verify.pl``; transcribed from
its ``replaceZeros``, ``interpolateTimes`` and ``errNorm``): both series
shifted by the column's offset and zeroed within its zerotol; the gold
interpolated linearly at the test's times (the gold must span them); the
RMS over time (trapezoids) of (gold - test) / (reltol |gold| + abstol),
zero where |gold - test| < absdifftol -- in units of reltol, so 1 or
below passes.  Defaults: reltol 0.01; abstol, absdifftol, zerotol 1e-12;
offset 0.  A netlist's ``*COMP name key=value ...`` line sets a column's,
matched by the column's exact name (`comp_tolerances`).

`crossing_times`: where a waveform crosses a level (mid-swing, for a
logic signal), linearly interpolated -- on a digital circuit the edges
dominate an RMS, and a delay says what moved.
"""
import re

import numpy as np

DEFAULTS = {'reltol': 0.01, 'abstol': 1e-12, 'absdifftol': 1e-12, 'zerotol': 1e-12,
            'offset': 0.0}


def read_prn(path):
    """A Xyce ``.prn``: (column names, rows as a float array)."""
    names, rows = None, []
    with open(path) as fh:
        for line in fh:
            s = line.strip()
            if not s:
                continue
            if names is None:
                ## (a column of an expression may hold blanks: `{I(V1) + 3E-3}`)
                names = re.findall(r'\{[^}]*\}|\S+', s)
            elif s.startswith('End of'):
                break
            else:
                rows.append([float(x) for x in s.split()])
    if names is None:
        raise ValueError(f'{path}: no header')
    return names, np.array(rows, dtype=float).reshape(len(rows), len(names))


def read_ngspice_print(text):
    """ngspice's batch ``.print`` output (the text): (column names, rows as
    a float array, one row per index, the columns in their first
    appearance's order)."""
    cols, values, current = [], {}, None
    for line in text.splitlines():
        words = line.split()
        if not words:
            continue
        if words[0] == 'Index' and len(words) >= 2:
            current = words[1:]
            for c in current:
                if c not in values:
                    cols.append(c)
                    values[c] = {}
            continue
        if current is None or not re.fullmatch(r'\d+', words[0]):
            continue
        if len(words) != len(current) + 1:
            continue
        k = int(words[0])
        for c, w in zip(current, words[1:], strict=True):
            values[c][k] = float(w)
    if not cols:
        raise ValueError('no ngspice .print table in the text')
    index = sorted(values[cols[0]])
    for c in cols:
        if sorted(values[c]) != index:
            raise ValueError(f'column {c} has other rows than {cols[0]}')
    return cols, np.array([[values[c][k] for c in cols] for k in index], dtype=float)


def read_ibm_output(lines):
    """The IBM power grid benchmarks' reference output (its lines):
    {node: (t, v)}."""
    out, node, rows = {}, None, []
    for line in lines:
        s = line.strip()
        if not s:
            continue
        if s.startswith('Node:'):
            if node is not None:
                out[node] = _pairs(rows)
            node, rows = s.split(None, 1)[1].strip(), []
        elif node is not None:
            rows.append([float(x) for x in s.split()])
    if node is not None:
        out[node] = _pairs(rows)
    return out


def _pairs(rows):
    a = np.array(rows, dtype=float).reshape(len(rows), 2)
    return a[:, 0], a[:, 1]


def comp_tolerances(comps):
    """The `*COMP` lines (their texts after ``*COMP``) as {column name:
    {key: value}} -- the keys ``reltol``, ``abstol``, ``absdifftol``,
    ``zerotol``, ``offset`` (``numfail``, a DC sweep's, is left)."""
    out = {}
    for text in comps:
        words = text.split()
        if not words:
            continue
        tol = out.setdefault(words[0], {})
        for w in words[1:]:
            k, _, v = w.partition('=')
            if k.lower() in DEFAULTS:
                tol[k.lower()] = float(v)
    return out


def _zeroed(v, offset, zerotol):
    v = np.asarray(v, dtype=float) + offset
    return np.where(np.abs(v) <= zerotol, 0.0, v)


def interpolate_at(t, gold_t, gold):
    """`gold` at the times `t` as ``interpolateTimes``: the first gold time
    at or after each `t` -- its value where the times are equal, else the
    line from the one before."""
    t, gold_t, gold = (np.asarray(a, dtype=float) for a in (t, gold_t, gold))
    if gold_t[0] > t[0] or gold_t[-1] < t[-1]:
        raise ValueError('the gold series does not span the test series')
    j = np.searchsorted(gold_t, t, side='left')
    lo = np.maximum(j - 1, 0)
    exact = gold_t[j] == t
    with np.errstate(divide='ignore', invalid='ignore'):
        line = (gold[j] - gold[lo]) / (gold_t[j] - gold_t[lo]) * (t - gold_t[lo]) + gold[lo]
    return np.where(exact, gold[j], line)


def xyce_verify(t, test, gold_t, gold, **tol):
    """The RMS relative error of one column, in units of its reltol (see
    the module note; 1 or below passes).  `tol`: the column's reltol,
    abstol, absdifftol, zerotol, offset (`DEFAULTS` for the rest)."""
    tol = {**DEFAULTS, **tol}
    unknown = set(tol) - set(DEFAULTS)
    if unknown:
        raise TypeError(f'unknown tolerances {sorted(unknown)}')
    t = np.asarray(t, dtype=float)
    f = _zeroed(test, tol['offset'], tol['zerotol'])
    g = interpolate_at(t, gold_t, _zeroed(gold, tol['offset'], tol['zerotol']))
    d = g - f
    scaled = d / (tol['reltol'] * np.abs(g) + tol['abstol'])
    if len(t) == 1:
        return float(abs(scaled[0]))
    integrand = np.where(np.abs(d) < tol['absdifftol'], 0.0, scaled * scaled)
    total = np.sum(0.5 * (integrand[1:] + integrand[:-1]) * np.abs(np.diff(t)))
    return float(np.sqrt(total / abs(t[-1] - t[0])))


def crossing_times(t, v, level, rising=None):
    """The times `v` crosses `level` (rising, falling, or either: None),
    linearly interpolated between the samples either side."""
    t, v = np.asarray(t, dtype=float), np.asarray(v, dtype=float)
    a, b = v[:-1] - level, v[1:] - level
    hit = (a < 0) & (b >= 0) if rising else (a > 0) & (b <= 0) if rising is False else (
        ((a < 0) & (b >= 0)) | ((a > 0) & (b <= 0)))
    k = np.flatnonzero(hit)
    return t[k] + (t[k + 1] - t[k]) * (-a[k]) / (b[k] - a[k])

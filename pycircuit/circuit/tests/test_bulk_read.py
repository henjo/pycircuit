"""THE BULK READ (`spice_import._Take`, plan 261 L21 B4): under `bulk` a
top-level line that is a bankable element on its own words alone is read
straight into its bank, with no `Card`, `Flat` or `Mapped` -- the circuit
the same, byte for byte, as reading every line as a record first
(`BULK_READ = False`): its nodes in order, its elements, each bank's arrays,
the merged nodes, the report; where that cannot be shown before building
(a name read twice, a pulse whose defaults need a missing `.tran`) every
line is read as a record.  `bulk=None` follows the file's size
(`BULK_FROM`); `write_ngspice` writes the banks from their arrays."""
from pathlib import Path

import pytest

from pycircuit.circuit import _bulk, spice_import
from pycircuit.utilities import spicecard

#: Every form the bulk read takes or leaves, its kind's banked records among
#: the lines it takes (the order), a subcircuit, an `.include`, a
#: continuation, a `.control` block and lines after `.end`.  `taken`: the
#: lines of each kind the read takes.
_DECK = """* the bulk read's forms
.tran 10p 2n
v1 1 0 1.8
r1 1 2 10
R2 2 N3 1K
r3 n3 4 r=20
c1 2 5 1p
C2 N3 0 2P
c3 4 0 1p ic=0.5
r4 4 5 {5*2}
r5 5 0 0
r6 (5 6) 30
i1 6 0 pulse(0 1m 0 10p 10p 0.5n 1n)
I2 N3 0 0 PULSE(0, 2m, 0.1n, 20p, 20p, 0.4n, 1n)
i3 4 0 dc 0 pulse(0 1m 0 10p 10p 0.5n 1n)
i4 n3 0 pulse(0 1m)
i5 2 0 0 ac 1 pulse(0 1m 0 10p 10p 0.5n 1n)
i6 6 0 1m
i7 6 0 pulse(0 1m 0 0 0 0 0)
r7 6 7
* a comment between
+ 1meg
r8 7 0 2
.subckt cell p q
rs p q 50
cs q 0 1p
.ends
x1 7 8 cell
r9 8 0 40
.include inc.cir
r11 9 0 70
.control
r99 1 0 1
.endc
{refused}.end
r98 1 0 5
"""

_INC = """r10 8 9 60
c4 9 0 3p
i8 9 0 pulse(0 1m 0 10p 10p 0.5n 1n)
"""

#: refused on either path (strict=False reports them), and a node that is an
#: expression (kept as written: the records' alone)
_REFUSED = """r.bad 1 0 10
r12 9 n.1 80
i9 9 0 pulse(0 1m 0 10p 10p 0.5n 1n) td=1n
i10 {N} 0 pulse(0 1m 0 10p 10p 0.5n 1n)
"""


def _write(tmp_path, refused=True, deck=None):
    p = tmp_path / 'deck.cir'
    (tmp_path / 'inc.cir').write_text(_INC)
    p.write_text((deck if deck is not None else _DECK).replace('{refused}', _REFUSED if refused else ''))
    return str(p)


def _snap(imp):
    """What an import builds: n, its nodes in order, its elements by name
    and class, each bank's kind, nodes and arrays (bytes), the merged
    nodes, the report, the .tran."""
    c = imp.circuit
    return {'n': c.n, 'nodes': [str(n.name) for n in c.nodes],
            'elements': [(k, type(e).__name__) for k, e in c.elements.items()],
            'banks': [(bk.name, bk.kind, [str(n.name) for n in bk.nodes],
                       *(None if x is None else (x.dtype.str, x.shape, x.tobytes())
                         for x in (bk.a, bk.b, bk.value, bk.wave)))
                      for bk in _bulk.banks(c)],
            'merged': dict(imp.merged), 'report': list(imp.report), 'tran': imp.tran}


def _both(monkeypatch, path, **kw):
    """The import read with records first (`BULK_READ = False`) and with
    the bulk read."""
    out = []
    for read in (False, True):
        monkeypatch.setattr(spice_import, 'BULK_READ', read)
        out.append(spice_import.import_netlist(path, bulk=True, **kw))
    return out


@pytest.mark.parametrize('merge', [False, True])
def test_the_bulk_read_builds_the_records_circuit(tmp_path, monkeypatch, merge):
    """Every form: the circuit read with records first, byte for byte --
    the lines taken R 7 (upper case and a continuation among them, one from
    the `.include`), C 3, I 4 (a DC value and commas, zero edges from
    `.tran`); a parameter, an IC, an expression, brackets, a `dc` word, a
    pulse of two values, a subcircuit's, a node written as an expression
    and the refused (a '.', a source's parameter) left to the records, the
    banked ones among the taken in file order; r = 0 a source (merged
    away with `merge_shorts`, a taken capacitor's end with it), the
    `.control` block's and the lines after `.end` read by neither."""
    rec, got = _both(monkeypatch, _write(tmp_path), strict=False, merge_shorts=merge)
    assert got.bulk_read == {'R': 7, 'C': 3, 'I': 4}
    assert rec.bulk_read == {}
    assert len(rec.elements) - len(got.elements) == 14
    assert [len(bk) for bk in _bulk.banks(got.circuit)] == [11, 5, 7]
    assert '{N}' in got.circuit.nodenames and '{n}' not in got.circuit.nodenames
    assert sorted(got.circuit.elements) == sorted(['i5', 'i6', 'v1'] + ([] if merge else ['r5']))
    assert sum('r.bad' in r or 'r12' in r or 'i9' in r for r in got.report) == 3
    assert any('11 resistors, 5 capacitors, 7 pulse current sources' in r for r in got.report)
    assert _snap(got) == _snap(rec)


def test_the_banks_keep_the_file_order(tmp_path, monkeypatch):
    """The banked records placed among the taken lines where the file has
    them: the resistor bank's values in the deck's order (r3, r4, r6 and
    the subcircuit's rs records), the source bank's pulses (i3, i4)."""
    _, got = _both(monkeypatch, _write(tmp_path, refused=False))
    R, _C, I = _bulk.banks(got.circuit)
    assert R.value.tolist() == [10.0, 1e3, 20.0, 10.0, 30.0, 1e6, 2.0, 50.0, 40.0, 60.0, 70.0]
    assert I.wave[:, 1].tolist() == [1e-3, 2e-3, 1e-3, 1e-3, 1e-3, 1e-3]
    ## (i4's two values and i7's zero edges: .tran's TSTEP and TSTOP)
    tstep, tstop = spicecard.number('10p'), spicecard.number('2n')
    assert I.wave[3].tolist() == I.wave[4].tolist() == [0.0, 1e-3, 0.0, tstep, tstep, tstop, tstop]


_TWICE = """* a name twice
.tran 10p 2n
v1 1 0 1
r1 1 2 10
{second}
c1 2 0 1p
.end
"""


@pytest.mark.parametrize('second', ['r1 2 0 20', 'r1 2 0 r=20', 'R1 2 0 20'])
def test_a_name_read_twice_reads_every_line_as_a_record(tmp_path, monkeypatch, second):
    """A name taken twice, or taken and a record's: the record path's
    renaming (`r1#1`, its note), every line read as a record, the report
    saying so last."""
    rec, got = _both(monkeypatch, _write(tmp_path, deck=_TWICE.replace('{second}', second)))
    assert got.bulk_read == {}
    assert got.report[-1].endswith('bulk: every line read as a record (an element name read twice)')
    a, b = _snap(rec), _snap(got)
    assert b.pop('report')[:-1] == a.pop('report')
    assert a == b
    assert any('renamed r1#1' in r for r in rec.report)


_NO_TRAN = """* a pulse's defaults without a .tran
v1 1 0 1
r1 1 0 10
i1 1 0 pulse(0 1m 0 0 10p 1n 2n)
.end
"""


def test_a_pulse_needing_a_missing_tran_reads_every_line_as_a_record(tmp_path, monkeypatch):
    """A taken pulse whose rise comes from `.tran` in a netlist without one:
    strict, the record path's error; else its report, the reason last."""
    p = _write(tmp_path, deck=_NO_TRAN)
    for read in (False, True):
        monkeypatch.setattr(spice_import, 'BULK_READ', read)
        with pytest.raises(spice_import.SpiceImportError, match="PULSE's defaults come from .tran"):
            spice_import.import_netlist(p, bulk=True)
    rec, got = _both(monkeypatch, p, strict=False)
    assert got.report[-1].endswith("(a pulse's defaults need a .tran, which this netlist lacks)")
    a, b = _snap(rec), _snap(got)
    assert b.pop('report')[:-1] == a.pop('report')
    assert a == b


#: 0 V sources: plain lines (vs1 and VS2 in a chain, vs3 to ground, a value
#: written 0.0) and, left to the records, a `dc` word (vs4: a 0 V source
#: all the same), an AC value (vs5), nonzero values (v1, vs6), a '.' in a
#: node (vs7: refused on either path)
_SHORTS = """* 0 V sources
.tran 10p 2n
v1 1 0 1.8
vs1 1 1b 0
VS2 1B 1c 0.0
vs3 2 0 0
r1 1c 2 10
r2 1c 3 20
vs4 3 3b dc 0
vs5 3b 4 0 ac 1
vs6 4 5 1
c1 5 0 1p
r3 5 0 1k
vs7 5 n.5 0
{extra}.end
"""


def _v_cards(imp):
    return sorted(c.name for c in imp.netlist.cards if c.kind == 'v')


@pytest.mark.parametrize('merge', [False, True])
def test_a_0v_source_is_read_into_the_merge(tmp_path, monkeypatch, merge):
    """With `merge_shorts` a plain 0 V source's line is taken with no
    record -- its nodes into the merge, beside the records' 0 V sources --
    the records' circuit byte for byte (the merged nodes, the report's
    count of the sources left out); without, every source a record (an
    element: nothing taken)."""
    rec, got = _both(monkeypatch, _write(tmp_path, deck=_SHORTS.replace('{extra}', '')), strict=False,
                     merge_shorts=merge)
    assert _snap(got) == _snap(rec)
    assert got.bulk_read == {'R': 3, 'C': 1, 'I': 0}
    assert sum('vs7' in r for r in got.report) == 1
    every = ['v1', 'vs1', 'vs2', 'vs3', 'vs4', 'vs5', 'vs6', 'vs7']
    assert _v_cards(rec) == every
    assert _v_cards(got) == (['v1', 'vs4', 'vs5', 'vs6', 'vs7'] if merge else every)
    if merge:
        assert got.merged == {'1b': '1', '1c': '1', '2': '0', '3b': '3'}
        assert any('merge_shorts: 4 0 V sources left out, 4 nodes merged' in r for r in got.report)
        assert sorted(got.circuit.elements) == ['v1', 'vs5', 'vs6']


@pytest.mark.parametrize('extra, why', [
    ('.print tran i(vs1)\n', 'a 0 V source whose current a .print reads'),
    ('vs1 4 0 dc 0\n', 'an element name read twice')])
def test_a_0v_source_the_merge_cannot_take_reads_every_line_as_a_record(tmp_path, monkeypatch, extra, why):
    """A taken 0 V source whose current a `.print` reads (the merge keeps
    it: an element), or whose name a record has too: every line read as a
    record, the report saying why last."""
    rec, got = _both(monkeypatch, _write(tmp_path, deck=_SHORTS.replace('{extra}', extra)), strict=False,
                     merge_shorts=True)
    assert got.bulk_read == {}
    assert got.report[-1].endswith(f'bulk: every line read as a record ({why})')
    a, b = _snap(rec), _snap(got)
    assert b.pop('report')[:-1] == a.pop('report')
    assert a == b
    assert ('vs1' in got.circuit.elements) == extra.startswith('.print')


def test_bulk_follows_the_size_of_the_netlist(tmp_path, monkeypatch):
    """`bulk=None`: banks from `BULK_FROM` lines of the file (its lines
    counted, the boundary included), none below; `bulk=True` and `False`
    whatever its size."""
    p = _write(tmp_path, refused=False)
    lines = Path(p).read_text().count('\n')
    for limit, banked in ((lines, True), (lines + 1, False)):
        monkeypatch.setattr(spice_import, 'BULK_FROM', limit)
        imp = spice_import.import_netlist(p)
        assert bool(_bulk.banks(imp.circuit)) is banked, limit
        assert (imp.bulk_read == {'R': 7, 'C': 3, 'I': 4}) is banked
    monkeypatch.setattr(spice_import, 'BULK_FROM', 10 ** 9)
    assert _bulk.banks(spice_import.import_netlist(p, bulk=True).circuit)
    monkeypatch.setattr(spice_import, 'BULK_FROM', 1)
    assert not _bulk.banks(spice_import.import_netlist(p, bulk=False).circuit)


def _elements(imp, nodes=None):
    """An element-built import's elements by class, nodes (through
    `nodes`, the written deck's) and parameters, sorted."""
    nodes = nodes or {}
    return sorted((m.cls.__name__, tuple(nodes.get(n, n) for n in m.nodes),
                   tuple(sorted(m.params.items()))) for m in imp.elements)


def test_write_ngspice_writes_the_banks_from_their_arrays(tmp_path, monkeypatch):
    """A bulk import written: its deck read back one element each gives the
    element-built import's elements -- classes, nodes, values bit for bit
    -- each once (a banked record is the bank's, not written twice); a
    banked element has no current to print."""
    p = _write(tmp_path, refused=False)
    imp = spice_import.import_netlist(p, bulk=True)
    assert imp.bulk_read == {'R': 7, 'C': 3, 'I': 4}
    out = str(tmp_path / 'written.cir')
    names = imp.write_ngspice(out)
    back = spice_import.import_netlist(out, bulk=False)
    ref = spice_import.import_netlist(p, bulk=False)
    assert _elements(back) == _elements(ref, names)
    text = Path(out).read_text()
    assert 'rbank10 ' in text and 'cbank4 ' in text and 'ibank5 ' in text
    with pytest.raises(spice_import.SpiceImportError, match='no element r1'):
        imp.write_ngspice(out, probes=[('i', 'r1')])


def test_a_bank_name_stays_apart_from_an_elements(tmp_path):
    """An element whose written name begins with a bank's prefix (a DC
    current source ``ibank0``): that bank's lengthened (``ibank_``), every
    name unique, each element read back."""
    p = _write(tmp_path, deck='* t\n.tran 10p 1n\nv1 1 0 1\nibank0 1 0 1m\n'
                              'i2 1 0 pulse(0 1m 0 10p 10p 0.5n 1n)\nr1 1 0 10\n.end\n')
    imp = spice_import.import_netlist(p, bulk=True)
    assert imp.bulk_read == {'R': 1, 'C': 0, 'I': 1}
    out = str(tmp_path / 'written.cir')
    imp.write_ngspice(out)
    words = [ln.split()[0] for ln in Path(out).read_text().splitlines() if ln and ln[0] in 'riv']
    assert sorted(words) == ['ibank0', 'ibank_0', 'rbank0', 'v1']
    back = spice_import.import_netlist(out, bulk=False)
    assert sorted(m.cls.__name__ for m in back.elements) == ['IPulse', 'IS', 'R', 'VS']

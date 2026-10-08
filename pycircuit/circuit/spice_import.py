"""Import a SPICE netlist: `import_netlist(path)` -> `Imported`.

The netlist is read by `pycircuit.utilities.spicenetlist` and flattened:
every subcircuit instance expanded in place, its elements and internal
nodes named by the instance path joined with ``:`` -- Xyce's separator,
so a gold column such as ``V(X1:3)`` names the node; never ``.``, which is
pycircuit's hierarchy separator.  SPICE's node ``0`` is ``gnd``.  Each
element is mapped as SPICE defines it:

====  ===============================================================
R     `R`; ``r = 0`` a 0 V source (SPICE's ammeter idiom)
C, L  `C`, `L`; an ``IC=`` only under UIC (SPICE ignores it otherwise)
K     `CoupledInductors`, replacing its two inductors
.ic   under UIC (or NOOP) the transient's starting values; without, the
      operating point solved with those nodes held, released at t = 0
.nodeset  the operating point solved from them (held, then released)
V, I  `VS`, `IS` and the waveform sources, SPICE's defaults filled in
      from ``.tran`` (PULSE's TR and TF a zero or absent one TSTEP, PW
      and PER TSTOP; SIN's FREQ 1/TSTOP; EXP's TAU TSTEP and TD2
      TD1 + TSTEP); the AC magnitude 0 unless ``AC`` is given; with a
      waveform the DC value is the waveform's (SPICE's transient)
G     `VCCS`, its control pins first
D     `DiodeSpiceHdl`
Q     `GummelPoonNpnHdl` / `GummelPoonPnpHdl`; a substrate node only
      where no substrate junction is asked for (``cjs = 0``)
M     `MosLevel1Hdl` / `MosLevel3Hdl` and their PMOS twins by LEVEL: a
      PMOS threshold as its magnitude, KP from UO and TOX where KP is
      not given, PHI 0.6 where neither PHI nor NSUB is given (as SPICE);
      TPG and NSS without effect where VTO is given, a VTO SPICE would
      derive from NSUB refused
X     flattened
====  ===============================================================

What cannot be mapped is refused: one `SpiceImportError` lists every
occurrence with its file and line (``strict=False``: into `report`
instead, the element left out).  What is read but has no effect here --
an ``OFF`` flag, an option group, a ``.print`` -- goes into `report`, as
does a second element of one flat name (renamed ``name#1``).  A name
with a ``.`` is refused.  The SPICE benchmark plan's stage 2.
"""
import re
from dataclasses import dataclass

from pycircuit.circuit import elements, elements_hdl
from pycircuit.circuit.circuit import SubCircuit, defaultepar, gnd
from pycircuit.utilities import spicecard, spicenetlist
from pycircuit.utilities.spicecard import SpiceCardError
from pycircuit.utilities.spicenetlist import SpiceNetlistError

#: The permittivity SPICE's level 1 takes for the oxide (mos1temp.c).
EPS_OX = 3.9 * 8.854214871e-12


class SpiceImportError(SpiceNetlistError):
    """What the importer cannot map -- every occurrence, with its file and
    line."""


#: SPICE's name of a device parameter -> the pycircuit class's.
DIODE_PARAMS = {'is': 'IS', 'rs': 'rs', 'n': 'n', 'tt': 'tt', 'cjo': 'cjo', 'cj0': 'cjo',
                'cj': 'cjo', 'vj': 'vj', 'pb': 'vj', 'm': 'm', 'mj': 'm', 'eg': 'eg',
                'xti': 'xti', 'fc': 'fc', 'bv': 'bv', 'ibv': 'ibv', 'kf': 'kf', 'af': 'af',
                'tnom': 'tnom'}
BJT_PARAMS = {'is': 'IS', 'bf': 'bf', 'nf': 'nf', 'vaf': 'vaf', 'va': 'vaf', 'ikf': 'ikf',
              'ik': 'ikf', 'ise': 'ise', 'ne': 'ne', 'br': 'br', 'nr': 'nr', 'var': 'var',
              'vb': 'var', 'ikr': 'ikr', 'isc': 'isc', 'nc': 'nc', 'rb': 'rb', 'rbm': 'rbm',
              're': 're', 'rc': 'rc', 'cje': 'cje', 'vje': 'vje', 'pe': 'vje', 'mje': 'mje',
              'me': 'mje', 'tf': 'tf', 'xtf': 'xtf', 'vtf': 'vtf', 'itf': 'itf', 'cjc': 'cjc',
              'vjc': 'vjc', 'pc': 'vjc', 'mjc': 'mjc', 'mc': 'mjc', 'xcjc': 'xcjc', 'tr': 'tr',
              'fc': 'fc', 'xtb': 'xtb', 'eg': 'eg', 'xti': 'xti', 'kf': 'kf', 'af': 'af',
              'tnom': 'tnom'}
#: A bipolar card's substrate junction: mapped only where it is off.
BJT_SUBSTRATE = ('cjs', 'ccs', 'vjs', 'ps', 'mjs', 'ms')
MOS1_PARAMS = {'vto': 'vto', 'vt0': 'vto', 'kp': 'kp', 'gamma': 'gamma', 'phi': 'phi',
               'lambda': 'lambd', 'tox': 'tox', 'nsub': 'nsub', 'ld': 'ld', 'cgso': 'cgso',
               'cgdo': 'cgdo', 'cgbo': 'cgbo', 'cbd': 'cbd', 'cbs': 'cbs', 'is': 'IS',
               'pb': 'pb', 'cj': 'cj', 'cjsw': 'cjsw', 'mj': 'mj', 'mjsw': 'mjsw', 'fc': 'fc',
               'js': 'js', 'rd': 'rd', 'rs': 'rs', 'rsh': 'rsh', 'kf': 'kf', 'af': 'af',
               'tnom': 'tnom'}
MOS3_PARAMS = {'vto': 'vto', 'vt0': 'vto', 'kp': 'kp', 'uo': 'u0', 'u0': 'u0', 'gamma': 'gamma',
               'phi': 'phi', 'tox': 'tox', 'nsub': 'nsub', 'xj': 'xj', 'nfs': 'nfs',
               'eta': 'eta', 'delta': 'delta', 'theta': 'theta', 'vmax': 'vmax',
               'kappa': 'kappa', 'ld': 'ld', 'wd': 'wd', 'cgso': 'cgso', 'cgdo': 'cgdo',
               'cgbo': 'cgbo', 'cbd': 'cbd', 'cbs': 'cbs', 'is': 'IS', 'pb': 'pb', 'cj': 'cj',
               'cjsw': 'cjsw', 'mj': 'mj', 'mjsw': 'mjsw', 'fc': 'fc', 'js': 'js', 'rd': 'rd',
               'rs': 'rs', 'rsh': 'rsh', 'kf': 'kf', 'af': 'af', 'tnom': 'tnom'}
MOS_INSTANCE = ('l', 'w', 'ad', 'as', 'pd', 'ps', 'nrd', 'nrs')
MOS_CLASSES = {(1, 'nmos'): elements_hdl.MosLevel1Hdl, (1, 'pmos'): elements_hdl.MosLevel1PmosHdl,
               (3, 'nmos'): elements_hdl.MosLevel3Hdl, (3, 'pmos'): elements_hdl.MosLevel3PmosHdl}

#: A source's waveform: its pycircuit classes (V, I) and its parameters' names.
WAVES = {'pulse': (elements.VPulse, elements.IPulse, ('1', '2', 'td', 'tr', 'tf', 'pw', 'per')),
         'sin': (elements.VSin, elements.ISin, ('o', 'a', 'freq', 'td', 'theta', 'phase')),
         'exp': (elements.VExp, elements.IExp, ('1', '2', 'td1', 'tau1', 'td2', 'tau2')),
         'pwl': (elements.VPWL, elements.IPWL, None)}
_REFUSED_WAVES = ('sffm', 'am', 'pat', 'data', 'sweep')
_SOURCE_WORDS = ('dc', 'ac') + tuple(WAVES) + _REFUSED_WAVES

_KINDS = {'e': 'a voltage-controlled voltage source', 'f': 'a current-controlled current source',
          'h': 'a current-controlled voltage source', 'b': 'a behavioural source',
          's': 'a voltage-controlled switch', 'w': 'a current-controlled switch',
          't': 'a transmission line', 'o': 'a lossy transmission line', 'j': 'a JFET',
          'z': 'a MESFET', 'u': 'a distributed RC line', 'y': 'a Xyce device',
          'n': 'a digital device', 'a': 'a code model', 'p': 'a port'}


@dataclass
class Flat:
    """One element of the flattened netlist: its flat name, its card, the
    subcircuit it is defined in (None: the top level), the instance's
    node map and parameter values."""
    name: str
    card: object
    sub: object
    prefix: str
    ports: dict
    overrides: dict

    def node(self, word):
        if word in self.ports:
            return self.ports[word]
        return '0' if word == '0' else self.prefix + word


@dataclass
class Mapped:
    """An element as mapped: its flat name, SPICE kind and flat nodes, the
    pycircuit class and its keyword arguments, and (for a device) the model
    it was given -- what `write_ngspice` writes back."""
    name: str
    kind: str
    nodes: list
    cls: type
    params: dict
    model: object
    where: object


class Imported:
    """An imported netlist (see the module note): `circuit`, the flat
    `SubCircuit`; `elements`, the `Mapped` elements in order; `tran`, the
    `.tran` (``tstep``, ``tstop``, ``tstart``, ``tmax``, ``start`` --
    ``'op'``, ``'uic'`` or ``'noop'``) or None; `ic` and `nodeset`, flat
    node -> volts; `temp`, degrees C or None; `options`, the Transient's
    keywords the options ask for; `report`, what was read without effect
    (and, with ``strict=False``, what was refused); `netlist`, the read
    `spicenetlist.Netlist`."""

    def __init__(self, netlist):
        self.netlist = netlist
        self.circuit = None
        self.elements = []
        self.tran = None
        self.ic = {}
        self.nodeset = {}
        self.temp = None
        self.options = {}
        self.integrator = None
        self.merged = {}
        self.report = []

    def epar(self):
        """The analysis' environment: the netlist's temperature."""
        ep = defaultepar.copy()
        if self.temp is not None:
            ep.T = 273.15 + self.temp
        return ep

    def transient(self, **overrides):
        """The `.tran` as pycircuit's transient: ``(Transient, keyword
        arguments for its solve)``.  The netlist's options, its temperature,
        UIC/NOOP (zeros, plus `.ic`), `.ic` without UIC (the operating point
        with those nodes held), `.nodeset` and TMAX become the Transient's
        keywords; `overrides` replace any of them.  TSTART is not a
        Transient's: the run starts at 0 (slice its result)."""
        from pycircuit.circuit.transient import Transient
        if self.tran is None:
            raise SpiceImportError('this netlist has no .tran')
        kw = dict(self.options)
        kw['epar'] = self.epar()
        if self.integrator is not None:
            kw['integrator'] = self.integrator()
        if self.tran['tmax']:
            kw['timestep_max'] = self.tran['tmax']
        if self.tran['start'] in ('uic', 'noop'):
            kw['uic'] = True
        elif self.nodeset:
            kw['nodeset'] = dict(self.nodeset)
        if self.ic:
            ## (under UIC the starting values; without, held while the
            ## operating point is solved -- SPICE's `.ic`)
            kw['ic'] = dict(self.ic)
        kw.update(overrides)
        return (Transient(self.circuit, **kw),
                {'tend': self.tran['tstop'], 'timestep': self.tran['tstep']})

    def write_ngspice(self, path, probes=None):
        """The flat circuit as an ngspice deck at `path` -- every element
        with its values in full precision and SPICE's spelling, each
        device's model as mapped, the `.tran` and its UIC, `.ic`, the
        temperature and the options mapped, a `.print tran` of `probes`
        (flat nodes, and ``('i', element)`` for a source's current; default
        every node).  Names are made ngspice's (an element's starts with its
        kind; ``:`` becomes ``_``): returns the map from each flat node to
        its name in the deck, and from each current probe to its column in
        ngspice's output (``name#branch``)."""
        nodes = _NameMap()
        lines, models = [f'* {self.netlist.title} -- flat, written by pycircuit', ''], {}
        names = _NameMap()

        def name(m, kind):
            return names(m.name if m.name.startswith(kind) and ':' not in m.name
                         else f'{kind}_{m.name}')

        written = {}
        for m in self.elements:
            n = [nodes(x) for x in m.nodes]
            p = m.params
            written[m.name] = name(m, {elements.R: 'r', elements.C: 'c', elements.L: 'l',
                                       elements.CoupledInductors: 'k', elements.VCCS: 'g'}.get(
                m.cls) or _SOURCE_CLASSES.get(m.cls) or _DEVICE_CLASSES[m.cls][0])
            if m.cls is elements.R:
                lines.append(f'{written[m.name]} {n[0]} {n[1]} {p["r"]!r}')
            elif m.cls in (elements.C, elements.L):
                key = 'c' if m.cls is elements.C else 'L'
                ic = f' ic={p["ic"]!r}' if 'ic' in p else ''
                lines.append(f'{written[m.name]} {n[0]} {n[1]} {p[key]!r}{ic}')
            elif m.cls is elements.CoupledInductors:
                k = written[m.name]
                lines += [f'l{k}_1 {n[0]} {n[1]} {p["L1"]!r}', f'l{k}_2 {n[2]} {n[3]} {p["L2"]!r}',
                          f'{k} l{k}_1 l{k}_2 {p["K"]!r}']
            elif m.cls is elements.VCCS:
                lines.append(f'{written[m.name]} {n[2]} {n[3]} {n[0]} {n[1]} {p["gm"]!r}')
            elif m.cls in _SOURCE_CLASSES:
                lines.append(f'{written[m.name]} {n[0]} {n[1]} ' + _source_spec(m.cls, p))
            elif m.cls in _DEVICE_CLASSES:
                kind, mtype, table, inst = _DEVICE_CLASSES[m.cls]
                card = {k: v for k, v in p.items() if k not in inst}
                if mtype == 'pmos' and 'vto' in card:
                    card['vto'] = -card['vto']
                spec = ' '.join(f'{table.get(k, k)}={v!r}' for k, v in sorted(card.items()))
                if kind == 'm':
                    level = 3 if m.cls in (elements_hdl.MosLevel3Hdl,
                                           elements_hdl.MosLevel3PmosHdl) else 1
                    spec = f'level={level} {spec}'
                key = (mtype, spec)
                if key not in models:
                    models[key] = f'{kind}model{len(models) + 1}'
                ip = ' '.join(f'{inst[k]}={v!r}' for k, v in p.items() if k in inst)
                lines.append(f'{written[m.name]} {" ".join(n)} {models[key]} {ip}'.rstrip())
            else:
                raise SpiceImportError(f'{m.name}: no ngspice form for {m.cls.__name__}')
        lines.append('')
        for (mtype, spec), mname in models.items():
            lines.append(f'.model {mname} {mtype} ({spec})')
        if self.tran is not None:
            t = self.tran
            uic = ' uic' if t['start'] in ('uic', 'noop') else ''
            lines.append(f'.tran {t["tstep"]!r} {t["tstop"]!r} 0 '
                         f'{t["tmax"] or min(t["tstep"], t["tstop"] / 50)!r}{uic}')
        if self.ic:
            lines.append('.ic ' + ' '.join(f'v({nodes(k)})={v!r}' for k, v in self.ic.items()))
        if self.temp is not None:
            lines.append(f'.temp {self.temp!r}')
        if 'reltol' in self.options:
            lines.append(f'.options reltol={self.options["reltol"]!r}')
        columns, shown = {}, []
        for pr in (probes if probes is not None else list(nodes.names)):
            if isinstance(pr, tuple):
                if pr[1] not in written:
                    raise SpiceImportError(f'no element {pr[1]} to print the current of')
                columns[pr] = f'{written[pr[1]]}#branch'
                shown.append(f'i({written[pr[1]]})')
            elif pr != '0':
                shown.append(f'v({nodes(pr)})')
        if self.tran is not None and shown:
            lines.append('.print tran')
            for k in range(0, len(shown), 8):
                lines.append('+ ' + ' '.join(shown[k:k + 8]))
        lines.append('.end')
        with open(path, 'w') as fh:
            fh.write('\n'.join(lines) + '\n')
        out = dict(nodes.names)
        out.update(columns)
        return out


class _NameMap:
    """Flat names -> ngspice's: ``:`` -> ``_``, made unique."""

    def __init__(self):
        self.names, self.taken = {}, set()

    def __call__(self, flat):
        if flat not in self.names:
            base = cand = flat.replace(':', '_')
            k = 1
            while cand in self.taken:
                cand, k = f'{base}_{k}', k + 1
            self.names[flat] = cand
            self.taken.add(cand)
        return self.names[flat]


#: The source classes, their SPICE kinds and waveforms' spellings.
_SOURCE_CLASSES = {elements.VS: 'v', elements.IS: 'i'}
for _v, _i, _names in WAVES.values():
    _SOURCE_CLASSES[_v], _SOURCE_CLASSES[_i] = 'v', 'i'


def _source_spec(cls, p):
    pre = 'v' if _SOURCE_CLASSES[cls] == 'v' else 'i'
    out = f'dc {p[pre]!r} ac {p[pre + "ac"]!r}'
    if 'phase' in p and cls in (elements.VS, elements.IS):
        out += f' {p["phase"]!r}'
    for wave, (v, i, names) in WAVES.items():
        if cls in (v, i):
            vals = p['tvpairs'] if names is None else [
                p[pre + k if k[0] in 'oa12' and len(k) == 1 else k] for k in names]
            out += f' {wave}(' + ' '.join(repr(x) for x in vals) + ')'
    return out


#: A device class: its SPICE kind and model type, its parameters' SPICE
#: names (where not its own), its instance parameters' SPICE names.
_DEVICE_CLASSES = {
    elements_hdl.DiodeSpiceHdl: ('d', 'd', {'IS': 'is'}, {'area': 'area'}),
    elements_hdl.GummelPoonNpnHdl: ('q', 'npn', {'IS': 'is'}, {'area': 'area'}),
    elements_hdl.GummelPoonPnpHdl: ('q', 'pnp', {'IS': 'is'}, {'area': 'area'}),
}
for _cls, _mtype in ((elements_hdl.MosLevel1Hdl, 'nmos'), (elements_hdl.MosLevel1PmosHdl, 'pmos'),
                     (elements_hdl.MosLevel3Hdl, 'nmos'), (elements_hdl.MosLevel3PmosHdl, 'pmos')):
    _DEVICE_CLASSES[_cls] = ('m', _mtype, {'IS': 'is', 'lambd': 'lambda', 'u0': 'uo'},
                             {k: k for k in MOS_INSTANCE} | {'asrc': 'as'})


def import_netlist(path, dialect='xyce', strict=True, merge_shorts=False):
    """Read and map the netlist `path` (see the module note).  `dialect`
    (``'xyce'`` or ``'ngspice'``) chooses where the two read a deck
    differently: the option and temperature statements.  `merge_shorts`:
    the nodes a 0 V source joins -- one whose current no `.print` reads --
    made one node and the source left out (a power grid's pads: ibmpg1t's
    unknowns roughly halve); `Imported.merged` maps each node merged away
    to the node it became."""
    if dialect not in ('xyce', 'ngspice'):
        raise ValueError(f'dialect {dialect!r}: xyce or ngspice')
    return _Importer(spicenetlist.read(path), dialect, strict, merge_shorts).run()


class _Importer:

    def __init__(self, net, dialect, strict, merge_shorts=False):
        self.net = net
        self.dialect = dialect
        self.strict = strict
        self.merge_shorts = merge_shorts
        self.out = Imported(net)
        self.errors = []
        self.inductors = {}         # flat name -> [Mapped, coupled by]
        self.couplings = []         # (Flat, flat L names, k)
        self.model_values = {}      # (model, instance parameters) -> values
        self.names = set()          # the flat element names taken

    ## -- reporting ----------------------------------------------------------

    def refuse(self, where, what):
        self.errors.append(f'{where}: {what}')

    def note(self, where, what):
        self.out.report.append(f'{where}: {what}')

    ## -- values -------------------------------------------------------------

    def value(self, raw, f):
        """A value of `f`'s card: a number, or an expression in the scope
        it is defined in, with the instance's parameters."""
        if spicecard._NUMBER.match(raw):
            return spicecard.number(raw)
        scope = f.sub.scope if f.sub is not None else self.net.deck.global_scope
        return self.net.deck.evaluate(raw, scope, **f.overrides)

    @staticmethod
    def is_value(word):
        return word[0] in "{'\"" or spicecard._NUMBER.match(word) is not None

    ## -- the run -------------------------------------------------------------

    def run(self):
        self.analyses()
        self.temperature()
        self.option_groups()
        for f in self.flatten():
            try:
                self.element(f)
            except SpiceCardError as e:
                self.refuse(f.card.where, f'{f.name}: {e}')
        self.couple()
        self.initial_conditions()
        if self.merge_shorts:
            self.merge()
        for what, where in self.net.unsupported:
            self.note(where, f'{what}: not read')
        for d in self.net.prints:
            self.note(d.where, f'.print {" ".join(d.words)}: not mapped (the probes are '
                               'the caller\'s)')
        if self.errors:
            if self.strict:
                raise SpiceImportError('this netlist cannot be imported:\n  '
                                       + '\n  '.join(self.errors))
            self.out.report = self.errors + self.out.report
        self.out.circuit = self.build()
        self.unholdable()
        return self.out

    def merge(self):
        """`merge_shorts` (see `import_netlist`): union-find over the 0 V
        sources no `.print` reads, ground the representative of its class;
        refused where a merge would short a source left in."""
        read = set()
        for d in self.net.prints:
            for w in d.words:
                read.update(m.lower() for m in re.findall(r'[iI]\(\s*([^,()\s]+)', w))
        parent = {}

        def find(n):
            while parent.get(n, n) != n:
                n = parent[n]
            return n

        shorts = [m for m in self.out.elements
                  if m.cls is elements.VS and m.params.get('v') == 0.0
                  and not m.params.get('vac') and not m.params.get('phase')
                  and m.name not in read]
        for m in shorts:
            a, b = find(m.nodes[0]), find(m.nodes[1])
            if a != b:
                keep, gone = (a, b) if a == '0' or (b != '0' and a < b) else (b, a)
                parent[gone] = keep
        dropped = {id(m) for m in shorts}
        self.out.elements = [m for m in self.out.elements if id(m) not in dropped]
        merged = {n: find(n) for n in parent}
        for m in self.out.elements:
            m.nodes = [merged.get(n, n) for n in m.nodes]
            if m.cls is elements.VS and m.nodes[0] == m.nodes[1]:
                self.refuse(m.where, f'{m.name}: merge_shorts would short this source')
        for target in (self.out.ic, self.out.nodeset):
            for n in list(target):
                if n in merged:
                    target.setdefault(merged[n], target.pop(n))
        self.out.merged = merged
        self.note(self.net.files[0], f'merge_shorts: {len(shorts)} 0 V sources left out, '
                                     f'{len(merged)} nodes merged')

    def unholdable(self):
        """An `.ic` (without UIC) on a node a voltage source or an inductor
        holds at DC has no effect in SPICE (the source wins): left out,
        and said."""
        tran = self.out.tran
        if not self.out.ic or tran is None or tran['start'] != 'op':
            return
        from pycircuit.circuit.dcanalysis import DC
        names = {str(n.name) for n in self.out.circuit.nodes}
        ## (with strict=False a refused element can take a node with it)
        ic = {n: v for n, v in self.out.ic.items() if n in names}
        if not ic or 'gnd' not in names:
            return
        for node in DC(self.out.circuit, epar=self.out.epar()).unholdable(ic):
            del self.out.ic[node]
            self.note(self.net.ic_where, f'.ic v({node}): a voltage source or an inductor '
                                         'holds the node at DC -- no effect, as in SPICE '
                                         '(the source wins)')

    def build(self):
        cir = SubCircuit()
        for m in self.out.elements:
            nodes = [gnd if n == '0' else n for n in m.nodes]
            cir[m.name] = m.cls(*nodes, **m.params)
        return cir

    def add(self, f, kind, nodes, cls, params, model=None):
        for n in [f.name] + nodes:
            if '.' in n:
                raise SpiceImportError(f'the name {n!r}: a "." is pycircuit\'s hierarchy '
                                       'separator')
        name, k = f.name, 1
        while name in self.names:
            name, k = f'{f.name}#{k}', k + 1
        if name != f.name:
            self.note(f.card.where, f'{f.name}: a second element of that name, renamed {name}')
        self.names.add(name)
        m = Mapped(name, kind, nodes, cls, params, model, f.card.where)
        self.out.elements.append(m)
        return m

    ## -- flattening -----------------------------------------------------------

    def flatten(self):
        """Every element of the netlist, subcircuit instances expanded."""
        out = []

        def walk(cards, sub, prefix, ports, overrides, depth):
            for card in cards:
                if card.kind != 'x':
                    out.append(Flat(prefix + card.name, card, sub, prefix, ports, overrides))
                    continue
                f = Flat(prefix + card.name, card, sub, prefix, ports, overrides)
                words = [w for w in card.words if w != 'params:']
                if not words:
                    self.refuse(card.where, f'{f.name}: names no subcircuit')
                    continue
                s = self.net.subckts.get(words[-1])
                if s is None:
                    self.refuse(card.where, f'{f.name}: no subcircuit {words[-1]}')
                    continue
                nodes = words[:-1]
                if len(nodes) != len(s.ports):
                    self.refuse(card.where, f'{f.name}: {len(nodes)} nodes for {s.name}\'s '
                                            f'{len(s.ports)} ports')
                    continue
                if depth > 64:
                    self.refuse(card.where, f'{f.name}: subcircuits nested past 64 (recursive?)')
                    continue
                try:
                    inst = {k: self.value(v, f) for k, v in card.params.items()}
                except SpiceCardError as e:
                    self.refuse(card.where, f'{f.name}: {e}')
                    continue
                walk(s.cards, s, f.name + ':', {p: f.node(n) for p, n in zip(s.ports, nodes)},
                     inst, depth + 1)

        walk(self.net.cards, None, '', {}, {}, 0)
        return out

    ## -- analyses, temperature, options -----------------------------------

    def analyses(self):
        trans = [d for d in self.net.analyses if d.name == 'tran']
        for d in self.net.analyses:
            if d.name not in ('tran', 'op'):
                self.note(d.where, f'.{d.name}: not mapped')
        if len(trans) > 1:
            self.refuse(trans[1].where, 'a second .tran')
        if not trans:
            return
        d = trans[0]
        words = list(d.words)
        start = 'op'
        if words and words[-1] in ('uic', 'noop'):
            start = words.pop()
        try:
            nums = [spicecard.number(w) if spicecard._NUMBER.match(w)
                    else self.net.deck.evaluate(w) for w in words]
        except SpiceCardError as e:
            self.refuse(d.where, f'.tran: {e}')
            return
        if not 2 <= len(nums) <= 4 or d.params:
            self.refuse(d.where, f'.tran {" ".join(d.words)}: TSTEP TSTOP [TSTART [TMAX]] '
                                 '[UIC|NOOP]')
            return
        nums += [0.0] * (4 - len(nums))
        self.out.tran = {'tstep': nums[0], 'tstop': nums[1], 'tstart': nums[2],
                         'tmax': nums[3], 'start': start}
        if nums[2]:
            self.note(d.where, f'.tran TSTART {nums[2]:g}: the run starts at 0 (slice it)')

    def temperature(self):
        temps = []
        for d in self.net.temps:
            temps += [(w, d.where) for w in d.words]
        for d in self.net.options:
            if ('temp' in d.params and (self.dialect == 'ngspice'
                                        or d.words[:1] == ['device'])):
                temps.append((d.params['temp'], d.where))
        if len(temps) > 1:
            self.refuse(temps[1][1], 'a second temperature (a sweep)')
        elif temps:
            self.out.temp = self.net.deck.evaluate(temps[0][0])

    def option_groups(self):
        """Xyce's ``.options timeint reltol= method=`` (ngspice's ungrouped
        ``.options reltol= method=``) become the Transient's; every other
        option is reported and left."""
        from pycircuit.circuit import integrator
        methods = {'gear': integrator.Gear2Integrator, 'bdf': integrator.Gear2Integrator,
                   'trap': integrator.TrapezoidalIntegrator,
                   'trapezoidal': integrator.TrapezoidalIntegrator,
                   'euler': integrator.EulerIntegrator, 'be': integrator.EulerIntegrator}
        for d in self.net.options:
            if self.dialect == 'xyce':
                group = d.words[0] if d.words else ''
                if group == 'device' and set(d.params) <= {'temp'}:
                    continue
                read = group == 'timeint'
            else:
                group, read = '', True
            for k, v in d.params.items():
                if read and k == 'reltol':
                    self.out.options['reltol'] = spicecard.number(v)
                elif read and k == 'method' and v.lower() in methods:
                    self.out.integrator = methods[v.lower()]
                elif not (self.dialect == 'ngspice' and k == 'temp'):
                    self.note(d.where, f'.options {group + " " if group else ""}{k}={v}: '
                                       'not mapped')

    def initial_conditions(self):
        tran = self.out.tran
        for target, src, what in ((self.out.ic, self.net.ic, '.ic'),
                                  (self.out.nodeset, self.net.nodeset, '.nodeset')):
            for node, raw in src.items():
                try:
                    target[node] = self.net.deck.evaluate(raw)
                except SpiceCardError as e:
                    self.refuse(getattr(self.net, what[1:] + '_where'), f'{what} v({node}): {e}')
        if '0' in self.out.ic or '0' in self.out.nodeset:
            self.refuse(self.net.ic_where or self.net.nodeset_where,
                        'an .ic or .nodeset on the ground node')
        if self.net.nodeset and tran is not None and tran['start'] != 'op':
            self.note(self.net.nodeset_where, '.nodeset ignored under UIC (no operating '
                                              'point is solved)')

    ## -- elements -------------------------------------------------------------

    def element(self, f):
        kind = f.card.kind
        handler = getattr(self, 'map_' + kind, None)
        if handler is None:
            what = _KINDS.get(kind, f'an element of kind {kind!r}')
            self.refuse(f.card.where, f'{f.name}: {what} is not supported')
            return
        handler(f)

    def two_nodes(self, f):
        if len(f.card.words) < 2:
            raise SpiceImportError('two nodes needed')
        return [f.node(w) for w in f.card.words[:2]]

    def single_value(self, f, key):
        """A two-terminal's value: positional or ``key=``; any other
        parameter refused."""
        c = f.card
        params = dict(c.params)
        rest = c.words[2:]
        if key in params and not rest:
            raw = params.pop(key)
        elif len(rest) == 1 and key not in params:
            raw = rest[0]
        else:
            raise SpiceImportError(f'one value expected, read {" ".join(rest)!r} '
                                   f'{sorted(params)}')
        return self.value(raw, f), params

    def map_r(self, f):
        nodes = self.two_nodes(f)
        r, params = self.single_value(f, 'r')
        if params:
            raise SpiceImportError(f'resistor parameters {sorted(params)} are not supported')
        if r == 0:
            self.note(f.card.where, f'{f.name}: r = 0, a 0 V source')
            self.add(f, 'r', nodes, elements.VS, {'v': 0.0, 'vac': 0.0})
        else:
            self.add(f, 'r', nodes, elements.R, {'r': r})

    def _ic(self, f, params, kw):
        ic = params.pop('ic', None)
        if ic is None:
            return
        if self.out.tran is not None and self.out.tran['start'] in ('uic', 'noop'):
            kw['ic'] = self.value(ic, f)
        else:
            self.note(f.card.where, f'{f.name}: IC= ignored without UIC (as SPICE)')

    def map_c(self, f):
        nodes = self.two_nodes(f)
        c, params = self.single_value(f, 'c')
        kw = {'c': c}
        self._ic(f, params, kw)
        if params:
            raise SpiceImportError(f'capacitor parameters {sorted(params)} are not supported')
        self.add(f, 'c', nodes, elements.C, kw)

    def map_l(self, f):
        nodes = self.two_nodes(f)
        ind, params = self.single_value(f, 'l')
        kw = {'L': ind}
        self._ic(f, params, kw)
        if params:
            raise SpiceImportError(f'inductor parameters {sorted(params)} are not supported')
        self.inductors[f.name] = [self.add(f, 'l', nodes, elements.L, kw), None]

    def map_k(self, f):
        c = f.card
        if len(c.words) != 3 or c.params:
            raise SpiceImportError('K L1 L2 coupling expected')
        self.couplings.append((f, [f.prefix + w for w in c.words[:2]],
                               self.value(c.words[2], f)))

    def couple(self):
        for f, names, k in self.couplings:
            pair = []
            for n in names:
                rec = self.inductors.get(n)
                if rec is None:
                    self.refuse(f.card.where, f'{f.name}: no inductor {n}')
                elif rec[1] is not None:
                    self.refuse(f.card.where, f'{f.name}: {n} is coupled by {rec[1]} already '
                                              '(a coupling of three inductors)')
                else:
                    pair.append(rec)
            if len(pair) != 2:
                continue
            (a, _), (b, _) = pair
            if 'ic' in a.params or 'ic' in b.params:
                self.refuse(f.card.where, f'{f.name}: a coupled inductor\'s IC= is not supported')
                continue
            for rec in pair:
                rec[1] = f.name
                self.out.elements.remove(rec[0])
                self.names.discard(rec[0].name)
            self.add(f, 'k', a.nodes + b.nodes, elements.CoupledInductors,
                     {'L1': a.params['L'], 'L2': b.params['L'], 'K': k})

    def map_v(self, f):
        self.source(f, 0)

    def map_i(self, f):
        self.source(f, 1)

    def source(self, f, which):
        c = f.card
        nodes = self.two_nodes(f)
        if c.params:
            raise SpiceImportError(f'source parameters {sorted(c.params)} are not supported')
        words, i = c.words[2:], 0
        dc, ac, acphase, wave, args = None, None, 0.0, None, []
        while i < len(words):
            w = words[i]
            if w == 'dc' and i + 1 < len(words) and self.is_value(words[i + 1]):
                dc, i = self.value(words[i + 1], f), i + 2
            elif w == 'ac':
                ac, i = 1.0, i + 1
                if i < len(words) and self.is_value(words[i]):
                    ac, i = self.value(words[i], f), i + 1
                    if i < len(words) and self.is_value(words[i]):
                        acphase, i = self.value(words[i], f), i + 1
            elif w in WAVES and wave is None:
                wave, i = w, i + 1
                while i < len(words) and self.is_value(words[i]):
                    args.append(self.value(words[i], f))
                    i += 1
            elif w in _REFUSED_WAVES:
                raise SpiceImportError(f'a {w.upper()} waveform is not supported')
            elif self.is_value(w) and dc is None and wave is None:
                dc, i = self.value(w, f), i + 1
            else:
                raise SpiceImportError(f'cannot read {w!r}')
        pre = 'vi'[which]
        kw = {pre: 0.0 if wave else (dc or 0.0), pre + 'ac': ac or 0.0}
        if acphase:
            kw['phase'] = acphase
        if wave is None:
            self.add(f, c.kind, nodes, (elements.VS, elements.IS)[which], kw)
            return
        classes = WAVES[wave]
        if wave == 'pwl':
            if len(args) < 2 or len(args) % 2:
                raise SpiceImportError('PWL needs time-value pairs')
            kw['tvpairs'] = args
        else:
            wkw = self.waveform(wave, args, pre)
            if wkw.get('phase', 0.0) != kw.get('phase', wkw.get('phase', 0.0)):
                ## (one parameter is both: pycircuit's SIN phase is its AC phase)
                raise SpiceImportError('a SIN phase and a different AC phase')
            kw.update(wkw)
        self.add(f, c.kind, nodes, classes[which], kw)

    def waveform(self, wave, args, pre):
        """`wave`'s parameters, SPICE's defaults filled in from `.tran`."""
        names = WAVES[wave][2]
        if not 2 <= len(args) <= len(names):
            raise SpiceImportError(f'{wave.upper()} takes 2 to {len(names)} values')
        a = args + [0.0] * (len(names) - len(args))
        tran = self.out.tran
        needs = {'pulse': len(args) < 7 or 0.0 in args[3:7], 'sin': len(args) < 3 or a[2] == 0,
                 'exp': len(args) < 6 or 0.0 in (a[3], a[5]) or a[4] == 0}[wave]
        if needs and tran is None:
            raise SpiceImportError(f'{wave.upper()}\'s defaults come from .tran, which this '
                                   'netlist lacks')
        if wave == 'pulse':
            v1, v2, td, tr, tf, pw, per = a
            out = {'td': td, 'tr': tr or tran['tstep'], 'tf': tf or tran['tstep'],
                   'pw': pw or tran['tstop'], 'per': per or tran['tstop']}
            out[pre + '1'], out[pre + '2'] = v1, v2
        elif wave == 'sin':
            vo, va, freq, td, theta, phase = a
            out = {'freq': freq or 1.0 / tran['tstop'], 'td': td, 'theta': theta, 'phase': phase}
            out[pre + 'o'], out[pre + 'a'] = vo, va
        else:
            v1, v2, td1, tau1, td2, tau2 = a
            step = tran['tstep'] if tran else 0.0
            out = {'td1': td1, 'tau1': tau1 or step, 'td2': td2 or td1 + step,
                   'tau2': tau2 or step}
            out[pre + '1'], out[pre + '2'] = v1, v2
        return out

    def map_g(self, f):
        c = f.card
        if len(c.words) != 5 or c.params:
            raise SpiceImportError('G n+ n- nc+ nc- gm expected (a VALUE or POLY form is '
                                   'not supported)')
        n = [f.node(w) for w in c.words[:4]]
        self.add(f, 'g', [n[2], n[3], n[0], n[1]], elements.VCCS,
                 {'gm': self.value(c.words[4], f)})

    ## -- devices ------------------------------------------------------------

    def model(self, f, name, types):
        m = self.net.model(name, f.sub)
        if m is None:
            raise SpiceImportError(f'no model {name}')
        if m.type not in types:
            raise SpiceImportError(f'model {name} is a {m.type}, not a {" or ".join(types)}')
        key = (id(m), tuple(sorted(f.overrides.items())))
        if key not in self.model_values:
            try:
                self.model_values[key] = self.net.deck.values(m, **f.overrides)
            except SpiceCardError as e:
                raise SpiceImportError(f'model {name}: {e}') from e
        return m, dict(self.model_values[key])

    @staticmethod
    def level(values, name):
        level = values.pop('level', 1.0)
        if level != int(level):
            raise SpiceImportError(f'model {name}: LEVEL {level:g}')
        return int(level)

    @staticmethod
    def translate(values, table, name, device):
        out, left = {}, {}
        for k, v in values.items():
            if k not in table:
                left[k] = v
            elif table[k] in out:
                raise SpiceImportError(f'model {name}: {table[k]} given twice (as {k} too)')
            else:
                out[table[k]] = v
        if left:
            raise SpiceImportError(f'model {name}: {", ".join(sorted(left))} not supported by '
                                   f'{device}')
        return out

    def instance(self, f, words, allowed):
        """An instance's positional area and its ``name=value`` parameters
        (`allowed`); ``off`` is reported (it steers only SPICE's operating
        point search), anything else refused."""
        c = f.card
        rest = [w for w in words if w != 'off']
        if len(rest) != len(words):
            self.note(c.where, f'{f.name}: OFF ignored (an operating point hint)')
        params = dict(c.params)
        if 'm' in params and self.value(params['m'], f) == 1:
            params.pop('m')
        out = {}
        for k in list(params):
            if k in allowed:
                out[k] = self.value(params.pop(k), f)
        if 'ic' in params:
            self.note(c.where, f'{f.name}: IC= ignored (device initial conditions are not '
                               'supported)')
            params.pop('ic')
        if params:
            raise SpiceImportError(f'instance parameters {sorted(params)} are not supported')
        if len(rest) > 1 or (rest and not self.is_value(rest[0])):
            raise SpiceImportError(f'cannot read {" ".join(rest)!r}')
        if rest:
            if 'area' in out:
                raise SpiceImportError('an area given twice')
            out['area'] = self.value(rest[0], f)
        return out

    def map_d(self, f):
        c = f.card
        if len(c.words) < 3:
            raise SpiceImportError('D n+ n- model expected')
        nodes = self.two_nodes(f)
        m, values = self.model(f, c.words[2], ('d',))
        if self.level(values, m.name) != 1:
            raise SpiceImportError(f'model {m.name}: only diode LEVEL 1 is supported')
        kw = self.translate(values, DIODE_PARAMS, m.name, 'DiodeSpiceHdl')
        kw.update(self.instance(f, c.words[3:], ('area',)))
        self.add(f, 'd', nodes, elements_hdl.DiodeSpiceHdl, kw, model=m)

    def map_q(self, f):
        c = f.card
        w = c.words
        nn = next((k for k in (3, 4) if len(w) > k and self.net.model(w[k], f.sub) is not None),
                  None)
        if nn is None:
            raise SpiceImportError('Q c b e [s] model expected (no model found)')
        m, values = self.model(f, w[nn], ('npn', 'pnp'))
        if self.level(values, m.name) != 1:
            raise SpiceImportError(f'model {m.name}: only bipolar LEVEL 1 (Gummel-Poon) is '
                                   'supported')
        sub = {k: values.pop(k) for k in BJT_SUBSTRATE if k in values}
        if sub.get('cjs', 0.0) or sub.get('ccs', 0.0):
            raise SpiceImportError(f'model {m.name}: a substrate junction (CJS) is not '
                                   'supported')
        for k in ('irb', 'ptf'):
            if values.get(k, 0.0):
                raise SpiceImportError(f'model {m.name}: {k.upper()} is not supported')
            values.pop(k, None)
        kw = self.translate(values, BJT_PARAMS, m.name, 'GummelPoonHdl')
        kw.update(self.instance(f, w[nn + 1:], ('area',)))
        cls = elements_hdl.GummelPoonNpnHdl if m.type == 'npn' else elements_hdl.GummelPoonPnpHdl
        if nn == 4:
            self.note(c.where, f'{f.name}: the substrate node {w[3]} is left unconnected '
                               '(no substrate junction)')
        self.add(f, 'q', [f.node(x) for x in w[:3]], cls, kw, model=m)

    def map_m(self, f):
        c = f.card
        w = c.words
        if len(w) < 5:
            raise SpiceImportError('M d g s b model expected')
        m, values = self.model(f, w[4], ('nmos', 'pmos'))
        level = self.level(values, m.name)
        cls = MOS_CLASSES.get((level, m.type))
        if cls is None:
            stage = ' (the plan\'s stage 9)' if level == 2 else ''
            raise SpiceImportError(f'model {m.name}: MOS LEVEL {level} is not supported{stage}')
        geometry = {k: values.pop(k) for k in ('l', 'w') if k in values}
        if level == 1:
            uo = values.pop('uo', values.pop('u0', 600.0))
            if 'kp' not in values and 'tox' in values:
                values['kp'] = uo * 1e-4 * EPS_OX / values['tox']
        ## (SPICE derives VTO from NSUB -- with TPG and NSS -- only where VTO
        ## is not given; these classes take VTO as given, default 0)
        if 'nsub' in values and 'vto' not in values and 'vt0' not in values:
            raise SpiceImportError(f'model {m.name}: VTO from NSUB is not supported (give VTO)')
        for k in ('tpg', 'nss'):
            values.pop(k, None)
        kw = self.translate(values, MOS1_PARAMS if level == 1 else MOS3_PARAMS, m.name,
                            cls.__name__)
        if 'vto' in kw and m.type == 'pmos':
            kw['vto'] = -kw['vto']
        if 'phi' not in kw and 'nsub' not in kw:
            kw['phi'] = 0.6
        inst = self.instance(f, w[5:], MOS_INSTANCE)
        if 'area' in inst:
            raise SpiceImportError('a MOSFET takes no area')
        for k in ('l', 'w'):
            if k not in inst and k in geometry:
                inst[k] = geometry[k]
        if level == 1 and 'as' in inst:
            inst['asrc'] = inst.pop('as')
        kw.update(inst)
        self.add(f, 'm', [f.node(x) for x in w[:4]], cls, kw, model=m)


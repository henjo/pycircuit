"""The circuits at size (roadmap 282).  `chain`: the pulsed MOS chain of the earlier large-circuit gates, on the MOS's SPICE
defaults, whose stage gain is below one: 4 of its 4 005 unknowns swing.  `achain`: the same with a gain above one, so a pulse
travels it.  `mesh`: rows x cols nodes, each an inverter driven by the node above and joined to its four neighbours by resistors:
a 2-D pattern, so a factor with fill."""
from pycircuit.circuit import elements_hdl as eh
from pycircuit.circuit.elements import VS, C, R, SubCircuit, VPulse, gnd

ACTIVE = {'vto': 0.5, 'kp': 2e-4}


def _head():
    c = SubCircuit()
    c.add_node('vdd')
    c['vdd'] = VS('vdd', gnd, v=1.8)
    c.add_node('g0')
    c['vin'] = VPulse('g0', gnd, v1=0.0, v2=1.8, td=5e-9, tr=2e-9, tf=2e-9, pw=4e-8, per=1e-7)
    return c


def chain(n, mos=None):
    c = _head()
    for k in range(n):
        c.add_node(f'd{k}')
        c[f'rl{k}'] = R('vdd', f'd{k}', r=2e4)
        c[f'M{k}'] = eh.MosLevel1Hdl(f'd{k}', 'g0' if k == 0 else f'd{k - 1}', gnd, gnd, **(mos or {}))
        c[f'c{k}'] = C(f'd{k}', gnd, c=2e-14)
    return c


def mesh(rows, cols, rc=3e5, mos=ACTIVE):
    c = _head()
    def nm(i, j):
        return f'n{i}_{j}'
    for i in range(rows):
        for j in range(cols):
            c.add_node(nm(i, j))
    for i in range(rows):
        for j in range(cols):
            n = nm(i, j)
            c[f'rl{i}_{j}'] = R('vdd', n, r=2e4)
            c[f'M{i}_{j}'] = eh.MosLevel1Hdl(n, 'g0' if i == 0 else nm(i - 1, j), gnd, gnd, **mos)
            c[f'c{i}_{j}'] = C(n, gnd, c=2e-14)
            if j + 1 < cols:
                c[f'rh{i}_{j}'] = R(n, nm(i, j + 1), r=rc)
            if i + 1 < rows:
                c[f'rv{i}_{j}'] = R(n, nm(i + 1, j), r=rc)
    return c


def make(kind, size):
    if kind == 'chain':
        return chain(int(size))
    if kind == 'achain':
        return chain(int(size), ACTIVE)
    if kind == 'mesh':
        r, c = (int(v) for v in str(size).split('x'))
        return mesh(r, c)
    raise SystemExit(f'no circuit {kind}')

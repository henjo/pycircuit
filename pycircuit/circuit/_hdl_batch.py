"""ONE C CALL PER CLASS PER PASS: the stamp plan's batches of C-bound hdl
elements (2026-10-03; speed round 6).

Rounds 2 and 5 made every device evaluation of a chained hdl class a C
kernel; what a mid-sized model's pass was made of after them was the CALL
around each kernel.  Measured on a 20-MosLevel1 chain: `cir.G` 56 us, of
which the kernels' C work was 6.7 -- the other 50 us were twenty wrapped
calls (`el.G` -> the generated method -> `_chained_eval` -> `CKernel.
__call__`, 2.0 us each) through the plan's per-element loop (0.8-1.1 us
each).  Here the plan calls ONE driver per class per pass: a loop in C over
the class's elements calling the kernel's own function pointer -- the very
machine code the per-element path calls, on the same bytes:

* `X = x[NM]` is every element's `x[nm]` stacked (a fresh contiguous copy,
  as `x[nm]` is); each kernel reads its row.
* `PR[e]` points into the element's OWN pack (`_hdl_cp`, `CKernel.pack`),
  and the driver writes the temperature into its slot exactly as
  `CKernel.__call__` does (`p[t_index] = T`) -- so a pack is left as the
  per-element path leaves it, and a pass through either path reads the
  same bytes.  The batch MIRRORS whatever tuple the element holds at the
  pass and never repacks an existing one: `state_restore` (the branch
  check's speculative Newton) puts an old `__dict__` back with no epoch
  move, and the per-element path reads that old pack too.
* `OUT[e]` is written into the element's slots of the pass's buffer, so the
  `bincount` sums the same values in the same order.

Measured before it was built (the prototype, this box): the 20-MosLevel1
`G` pass 59.4 -> 7.5 us, `i` 64.9 -> 12.1, the Gummel-Poon chain the same,
20-PSP `G` -14 % and `i` -34 % (its C work dominates); every pass bit for
bit today's.

WHICH ELEMENTS, decided when the plan is built (`split`): a chained class
with no DC pins (the three that have them keep the per-element path's
`G_dc`/`i_dc` choice), the method the metaclass generated -- CODE identity
against the metaclass's nested definitions (`generated_code`), since a
spy, a `functools.wraps` spy and a Python override each carry their own
code object -- the circuit's own toolkit, a bound kernel whose `nx` is
the element's size and whose shape fills the element's slots; two or more
such elements of one class.  Batches are grouped by class, so a plain
subclass (it shares its parent's info and kernels) has a batch of its
own.

CHECKED ON EVERY PASS (`Batch.run`): `ENABLED`, the class still C-bound
and its kernel the one the batch cast (re-cast when it moved, refused
when its layout did), the method still the generated one (a class patch
during a run), the temperature one number (`CKernel`'s rule), and per
element: no instance shadow of the method (PCNR installs `el.i`/`el.G` on
the participating devices DURING a solve and removes them after) and a
pack that could be built.  An element that fails a per-element check is
handed back to the caller for the per-element path and the rest of its
class stays one call; a class that fails a class check takes the
per-element path for the pass.

THE SOURCE PASSES (`zero_source`, commit 2 of the round): every chained
library class's `u`, `u_dc` and `dudt` return a list of literal zeros (22
of 22 on 2026-10-03), and `cir.u` called each of them -- 44 calls at 1.3
us on the 20-MosLevel1 chain, 10 % of a step -- to add exact zeros to
bins that start at +0.0, which changes no bit.  `_add_element_subvectors`
skips such an element as it skips `Circuit.u`'s default, on the same
conditions (a numeric toolkit's scatter path, no instance shadow) plus no
`dtype` and not the 'ac' analysis (the complex `uac` path), decided once
per class by parsing the compiled function's source (`info['_u_zero']`,
a runtime key of the compile cache) and checked per call against the
method's code identity (a subclass override has its own code).

`ENABLED` (env `PYCIRCUIT_HDL_BATCH=0`) is read on every pass, so a
measurement can A/B it; `SKIP_ZERO_SOURCE` (`PYCIRCUIT_HDL_ZERO_U=0`) the
source passes' skip.  The driver is one small C source of its own,
exporting the chain functions' entry name with its own signature, built
and loaded once per process through `_hdl_cbackend.load_kernel` under its
own key (the limiter's precedent, `_hdl_climit`); a build that fails turns
batching off quietly (`STATUS`).

History: `doc/transient_history.md`, `_stamp_plan`; `doc/hdl_roadmap_260824.md`
sec. 63.
"""
import ast
import os
import types

import numpy as np

ENABLED = os.environ.get('PYCIRCUIT_HDL_BATCH', '1') != '0'
SKIP_ZERO_SOURCE = os.environ.get('PYCIRCUIT_HDL_ZERO_U', '1') != '0'

#: The pass driver.  `fn` is a chain kernel's function pointer (the object
#: `_dlopen` resolved, cast to `void *`); `X` is `n` rows of `sx` doubles,
#: `OUT` `n` rows of `so`; `PR[e]` the element's pack; `ti` the temperature
#: slot or -1 when the class reads none.
PASS_C = r"""
typedef void (*hdl_batch_fn_t)(const double *, const double *, double *);
void hdl_fn(const void *fn, const double *X, double **PR, double *OUT,
            long n, long sx, long so, long ti, double T)
{
    hdl_batch_fn_t f = (hdl_batch_fn_t) fn;
    long e;
    if (ti >= 0) {
        for (e = 0; e < n; e++) {
            PR[e][ti] = T;
            f(X + e * sx, PR[e], OUT + e * so);
        }
    } else {
        for (e = 0; e < n; e++)
            f(X + e * sx, PR[e], OUT + e * so);
    }
}
"""
PASS_CDEF = ('void hdl_fn(const void *fn, const double *X, double **PR, '
             'double *OUT, long n, long sx, long so, long ti, double T);')

#: `(ffi, cfn, double*, double**)` once the driver is loaded; False after a
#: build that failed (`STATUS` says why); None before the first use.
_driver = None
STATUS = 'not loaded'


def driver():
    """The loaded pass driver, or None where it cannot be had."""
    global _driver, STATUS
    if _driver is None:
        from pycircuit.circuit import _hdl_cbackend as cb
        try:
            ffi, cfn, _key, _cold, _secs = cb.load_kernel(PASS_C, PASS_CDEF)
        except (cb.CompileError, OSError) as e:
            _driver = False
            STATUS = f'off ({e})'
        else:
            _driver = (ffi, cfn, ffi.typeof('double *'), ffi.typeof('double **'))
            STATUS = 'c'
    return _driver or None


_GEN = None


def generated_code():
    """The code objects of the methods `BehaviouralMeta.__init__` defines
    for every class (`i`, `G`, `q`, `C`, `limit`, ...), by name -- each is
    defined exactly once there, so a class's method IS the generated one
    iff its `__code__` is that object."""
    global _GEN
    if _GEN is None:
        from pycircuit.circuit.hdl import BehaviouralMeta
        _GEN = {c.co_name: c for c in BehaviouralMeta.__init__.__code__.co_consts
                if isinstance(c, types.CodeType)}
    return _GEN


def is_generated(cls, m):
    """Whether `cls`'s method `m` is the one the metaclass generated."""
    code = getattr(getattr(cls, m, None), '__code__', None)
    return code is not None and code is generated_code().get(m)


class Batch:
    """One class's elements in one plan: the stacked node maps, the slots,
    the mirrored packs and the kernel's pointer (see the module note)."""

    __slots__ = ('NM', 'PR', 'cls', 'dst', 'els', 'entries', 'fnptr', 'info',
                 'k', 'kern', 'm', 'mirror', 'n', 'shape', 'slots', 'so', 'ti')

    def __init__(self, cls, m, entries, kern):
        self.cls, self.info, self.m = cls, cls._hdl_info, m
        self.entries = entries
        self.els = [e[1] for e in entries]
        self.NM = np.ascontiguousarray(
            np.stack([np.asarray(e[2], dtype=np.int64) for e in entries]))
        self.slots = [(e[3], e[4]) for e in entries]
        self.dst = np.concatenate([np.arange(a, b) for a, b in self.slots])
        self.n = len(entries)
        self.k, self.shape = kern.nx, kern.shape
        self.so = int(np.prod(kern.shape))
        self.mirror = [None] * self.n
        self.PR = np.zeros(self.n, dtype=np.uintp)
        self.kern = self.fnptr = None
        self.ti = -1
        self._take(kern)

    def _take(self, kern):
        """Cast `kern`'s function pointer for the driver; False where the
        kernel is not one the batch was built for."""
        from pycircuit.circuit import _hdl_cbackend as cb
        drv = driver()
        if (drv is None or not isinstance(kern, cb.CKernel)
                or kern.nx != self.k or kern.shape != self.shape):
            self.kern = None
            return False
        self.kern = kern
        self.fnptr = drv[0].cast('void *', kern.cfn)
        self.ti = -1 if kern.t_index is None else kern.t_index
        return True

    def run(self, x, epar):
        """The kernel's outputs for this pass: `(OUT, None)` with a row per
        element, `(OUT, pos)` with a row per element in `pos` (the others
        are the caller's), or None where the whole class is the caller's."""
        if not ENABLED:
            return None
        drv = _driver
        if not drv:
            return None
        info = self.info
        if not info.get('_c_bound'):
            return None
        kern = info['funcs'][self.m].__dict__.get('_hdl_c')
        if kern is not self.kern and not self._take(kern):
            return None
        if not is_generated(self.cls, self.m):
            return None
        ti = self.ti
        T = 0.0
        if ti >= 0:
            T = getattr(epar, 'T', 300.0)
            if type(T) is not float:
                ## (an int or a 0-d array is one number -- `CKernel`'s rule;
                ## `p[t_index] = T` converts it as `float` does)
                if type(T) is not int and np.ndim(T) != 0:
                    return None
                try:
                    T = float(T)
                except (TypeError, ValueError, OverflowError):
                    return None
        els, mirror, PR, m = self.els, self.mirror, self.PR, self.m
        skip = None
        for e in range(self.n):
            d = els[e].__dict__
            if m in d:
                skip = (skip or []) + [e]
                continue
            cp = d.get('_hdl_cp')
            if cp is None:
                try:
                    cp = kern.pack(els[e])
                except (TypeError, ValueError):
                    cp = False
                d['_hdl_cp'] = cp
            if cp is False:
                skip = (skip or []) + [e]
                continue
            if cp is not mirror[e]:
                mirror[e] = cp
                PR[e] = cp[0].ctypes.data
        ffi, cfn, dptr, pptr = drv
        if skip is None:
            X = x[self.NM]
            OUT = np.empty((self.n,) + self.shape)
            cfn(self.fnptr, ffi.from_buffer(dptr, X), ffi.from_buffer(pptr, PR),
                ffi.from_buffer(dptr, OUT), self.n, self.k, self.so, ti, T)
            return OUT, None
        pos = [e for e in range(self.n) if e not in skip]
        if not pos:
            return None
        X = x[self.NM[pos]]
        PRp = np.ascontiguousarray(PR[pos])
        OUT = np.empty((len(pos),) + self.shape)
        cfn(self.fnptr, ffi.from_buffer(dptr, X), ffi.from_buffer(pptr, PRp),
            ffi.from_buffer(dptr, OUT), len(pos), self.k, self.so, ti, T)
        return OUT, pos


def _returns_zeros(fn):
    """Whether the chain-compiled `fn`'s source ends in a `return` of a
    list of literal zeros (`0`, `0.0`, their negatives: adding any of them
    to a bin that starts at +0.0 changes no bit)."""
    src = getattr(fn, '_src', None)
    if src is None:
        return False
    try:
        tree = ast.parse(src)
    except SyntaxError:
        return False
    if len(tree.body) != 1 or not isinstance(tree.body[0], ast.FunctionDef):
        return False
    ret = tree.body[0].body[-1]
    if not isinstance(ret, ast.Return) or not isinstance(ret.value, ast.List):
        return False
    for e in ret.value.elts:
        if isinstance(e, ast.UnaryOp) and isinstance(e.op, ast.USub):
            e = e.operand
        if not (isinstance(e, ast.Constant) and type(e.value) in (int, float)
                and e.value == 0):
            return False
    return True


def zero_source(cls, m):
    """Whether `cls`'s generated `m` ('u' or 'dudt') adds nothing on the
    scatter path: the class's compiled `u` and `u_dc` (or `dudt`) return
    literal zeros -- decided once per class (`info['_u_zero']`) -- and the
    method is the generated one (an override has its own code; the
    generated method's other answers are zeros too)."""
    info = getattr(cls, '_hdl_info', None)
    if info is None:
        return False
    z = info.get('_u_zero')
    if z is None:
        z = info['_u_zero'] = {}
    r = z.get(m)
    if r is None:
        funcs = info['funcs']
        names = ('u', 'u_dc') if m == 'u' else ('dudt',)
        r = z[m] = all(_returns_zeros(funcs.get(n)) for n in names)
    return r and is_generated(cls, m)


def split(cir, m, calls):
    """The plan's non-constant `calls` (`(inst, el, nm, a, b)` in dict
    order) split into `(rest, batches)`: the batches of the C-bound classes
    with two or more eligible elements, and the rest in their order."""
    from pycircuit.circuit import _hdl_cbackend as cb
    mat = m in ('G', 'C')
    groups = {}
    rest = []
    for c in calls:
        _inst, el, nm, a, b = c
        cls = type(el)
        info = getattr(cls, '_hdl_info', None)
        kern = None
        if (a >= 0 and info is not None and info.get('chained')
                and not info['state_meta']['dc_pins']
                and el.toolkit is cir.toolkit and m not in el.__dict__
                and is_generated(cls, m)):
            fn = info['funcs'].get(m)
            kern = fn.__dict__.get('_hdl_c') if fn is not None else None
            k = len(nm)
            if not (isinstance(kern, cb.CKernel) and kern.nx == k
                    and kern.shape == ((k, k) if mat else (k,))
                    and b - a == (k * k if mat else k)):
                kern = None
        if kern is None:
            rest.append(c)
        else:
            groups.setdefault(cls, (kern, []))[1].append(c)
    batches = []
    for cls, (kern, ents) in groups.items():
        if len(ents) >= 2:
            batches.append(Batch(cls, m, ents, kern))
        else:
            rest.extend(ents)
    if batches and len(rest) > 1:
        order = {id(c): j for j, c in enumerate(calls)}
        rest.sort(key=lambda c: order[id(c)])
    return rest, batches

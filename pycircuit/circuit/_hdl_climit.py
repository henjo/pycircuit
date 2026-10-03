"""THE hdl LIMITER IN C (speed round 5, 2026-10-03): a chained model's
generated `limit()` -- the probes' parameter chains, the laws, the
canonical orders and the device write-back -- as ONE C kernel per class,
the Python closure's answer bit for bit.

WHY.  After speed round 3 trimmed the closure (`hdl.py`, `BehaviouralMeta`:
the parameter cache, the lists, the reused ranking) `cir.limit` was still
33 % of a 20-MosLevel1 gear step: 18.8 us a call, of which the Python body
17.6 (lists, the two rankings, the sorts, the write-back), the
solution-reading parameter chain `von` 6.2, the laws 2.6.  A C kernel call
costs 1.7 us.  The parameter chains alone in C would recover a quarter of
that; the whole closure recovers nearly all of it.

WHAT IS PRINTED.  Every parameter function of the spec keeps its
ingredients (`_hdl_limit_par`: the chain statements it reaches, its
expression, its signature) -- the sympy form the C printer prints `i`,
`q`, `G`, `C` from (`hdl._CChainPrinter`, through the same symbol map,
`hdl._c_symmap`).  The statements of every probe's parameters are merged
into one list (one expression per symbol, so a symbol is computed once
and serves every probe with the same bits), the laws of `_limiting.py`
are transliterated below in Python's own forms (`max(a, b)` is `(b > a)
? b : a`, `abs` is `fabs`, `log` is libm's -- the backend's `-fno-builtin`
and `-ffp-contract=off` flags keep every operation where Python put it),
and the closure's algorithm follows statement by statement: the eager
parameter values at `x0`, `drift`, each group ranked by `(-|vlim - vn|,
ra, rb)` then applied with the shifts, the device write-back's Kruskal
and its anchors, the singles ranked the same way and written one end at a
time.  `(vlim - vin) * -1.0` is a multiplication, never a negation, and
`vin = (vorig + sa) - sb` is always computed: a NaN's sign bit and a
zero's sign follow Python's.

THE ONE THING C DOES NOT ANSWER.  Python's `sorted` and `min` on keys
holding a NaN order by timsort's and the set's accidents, which no other
sort reproduces (measured: a stable insertion sort disagrees in 29 % of
random NaN-keyed lists).  So the kernel DECLINES -- returns 1 and leaves
`out` as it came -- whenever a key it is about to sort or minimise is NaN
(a NaN or inf row of the state at a probe, a NaN parameter, infinities
meeting in the shifts), and the Python closure answers that call with its
own order and its own warnings.  With finite keys the orders are total
(the probes' row pairs are unique), and any stable sort is `sorted`.

THE CALL.  `int hdl_fn(const double *x, const double *p, double *out)`: the
chain functions' entry name and pointer layout, so the object is built,
stored and loaded by `_hdl_cbackend.load_kernel` like any other (keyed by
its own source, this prelude included); `x` is `x0`, the last accepted
point the parameter chains read, `out` ARRIVES as the copy of the
unlimited state and is limited in place.  `CLimitKernel` packs the
parameters as `CKernel` does (the element's shared `_hdl_cp`, dropped by
`update()`), serves a float, int or numpy-float temperature, and hands
anything unusual -- a short or non-numeric vector, a 0-d array, a
declined call -- back to the closure as `None`.

WHAT IT CHANGES.  Nothing in the answer.  A class whose limiter cannot be
printed (an unsupported node, `tanh` -- the backend's named ulp
exception -- a parameter without ingredients, a signature that is not
the packed layout) keeps its Python closure and its own backend status;
`_c_limit_status` says why.  The Python closure's parameter functions are
not called on a C-bound class -- a test counting them switches the kernel
off (`ENABLED`, `PYCIRCUIT_HDL_CLIMIT=0`) -- and a spec list mutated in place
AFTER the kernel was bound is not seen by it, as the stamp plan's contract
has it.

History: `doc/pss_log_260902.md`, 2026-10-03 (speed round 5).
"""
import os
import types

import numpy as np

from pycircuit.circuit import _hdl_cbackend as _cb

#: the law of each probe kind, as the C `_lim` dispatches it
KINDS = {'pnj': 0, 'fet': 1, 'vds': 2, 'delta': 3, 'id': 4}

#: THE SWITCH: False keeps the Python closure on every class (the chain
#: functions stay on C) -- `PYCIRCUIT_HDL_CLIMIT=0`, read at import; the
#: closure reads the attribute per call, so a test can flip it.  The
#: tests that count the closure's parameter functions do: a numpy pin of
#: the WHOLE class moved their recorded transients, since the 'auto'
#: Newton options read a numpy class differently (2026-10-03).
ENABLED = os.environ.get('PYCIRCUIT_HDL_CLIMIT', '1') != '0'


class Refused(Exception):
    """This class's limiter is not printed to C; the message says why."""


def limit_cdef():
    """The limiter's entry: the chain functions' name, an `int` return."""
    from pycircuit.circuit import hdl
    return f'int {hdl._C_ENTRY}(const double *x, const double *p, double *out);'


## ---------------------------------------------------------------------------
## The limiter's own C prelude: the laws of `_limiting.py` in Python's forms,
## Python's tuple order on the ranking keys, the device write-back.  Part of
## every limit kernel's source (so `source_key` covers it), never of
## `hdl._KERNEL_C` (which would rebuild every chain function's object).
## ---------------------------------------------------------------------------
_LIMIT_C = r"""
/* -- the hdl limiter: _limiting.py's laws in Python's own forms ----------- */
static double _pymax(double a, double b) { return (b > a) ? b : a; }
static double _pymin(double a, double b) { return (b < a) ? b : a; }

static double _lim_pnj(double vnew, double vold, double VT, double IS) {
  if (IS <= 0.0) return vnew;
  double vc = VT * log(VT / (IS * 1.414213562));
  if (fabs(vnew - vold) <= 2.0 * VT) return vnew;
  if (vnew > vc && vnew > 0.0) {
    if (vold > 0.0) {
      double arg = 1.0 + (vnew - vold) / VT;
      return (arg > 0.0) ? vold + VT * log(arg) : vc;
    }
    return VT * log(vnew / VT);
  }
  return vnew;
}

static double _lim_fet(double vnew, double vold, double vto) {
  double vtsthi = fabs(2.0 * (vold - vto)) + 2.0;
  double vtstlo = fabs(vold - vto) + 1.0;
  double vtox = vto + 3.5;
  double delv = vnew - vold;
  if (vold >= vto) {
    if (vold >= vtox) {
      if (delv <= 0.0) {
        if (vnew >= vtox) {
          if (-delv > vtstlo) vnew = vold - vtstlo;
        } else {
          vnew = _pymax(vnew, vto + 2.0);
        }
      } else {
        if (delv >= vtsthi) vnew = vold + vtsthi;
      }
    } else {
      if (delv <= 0.0) vnew = _pymax(vnew, vto - 0.5);
      else vnew = _pymin(vnew, vto + 4.0);
    }
  } else {
    if (delv <= 0.0) {
      if (-delv > vtsthi) vnew = vold - vtsthi;
    } else {
      double vtemp = vto + 0.5;
      if (vnew <= vtemp) {
        if (delv > vtstlo) vnew = vold + vtstlo;
      } else {
        vnew = vtemp;
      }
    }
  }
  return vnew;
}

static double _lim_vds(double vnew, double vold) {
  if (vold < 0.0) return -_lim_vds(-vnew, -vold);
  if (vold >= 3.5) {
    if (vnew > vold) vnew = _pymin(vnew, 3.0 * vold + 2.0);
    else if (vnew < 3.5) vnew = _pymax(vnew, 2.0);
  } else {
    if (vnew > vold) vnew = _pymin(vnew, 4.0);
    else vnew = _pymax(vnew, -0.5);
  }
  return vnew;
}

static double _lim_delta(double vnew, double vold, double vmax) {
  double d = vnew - vold;
  if (d > vmax) return vold + vmax;
  if (d < -vmax) return vold - vmax;
  return vnew;
}

/* `apply_limit`: pnj's pars are (IS, VT), the law takes (VT, IS) */
static double _lim(int kind, double vnew, double vold, double p0, double p1) {
  switch (kind) {
    case 0: return _lim_pnj(vnew, vold, p1, p0);
    case 1: return _lim_fet(vnew, vold, p0);
    case 2: return _lim_vds(vnew, vold);
    case 3: return _lim_delta(vnew, vold, p0);
    default: return vnew;
  }
}

/* Python's tuple `<` on (double, int, int) keys: the first element that
   is not `==` decides */
static int _tless(double a0, int a1, int a2, double b0, int b1, int b2) {
  if (!(a0 == b0)) return a0 < b0;
  if (a1 != b1) return a1 < b1;
  return a2 < b2;
}

/* `sorted(order, key=(k0[j], k1[j], k2[j]))` for finite k0: a stable
   insertion sort (the keys are indexed by the entries of `order`) */
static void _sort_by(int *order, int n, const double *k0, const int *k1,
                     const int *k2) {
  for (int i = 1; i < n; i++) {
    int key = order[i];
    int j = i - 1;
    while (j >= 0 && _tless(k0[key], k1[key], k2[key],
                            k0[order[j]], k1[order[j]], k2[order[j]])) {
      order[j + 1] = order[j];
      j--;
    }
    order[j + 1] = key;
  }
}

static int _find(int *parent, int n) {
  while (parent[n] != n) { parent[n] = parent[parent[n]]; n = parent[n]; }
  return n;
}

/* `device_writeback(out, targets, drift, pinned=moved)`: the rows it
   wrote are marked in `moved` AFTER the walk (the walk reads the set as it
   was).  Returns -1 where a Kruskal key is NaN (declined). */
#define _WB_MAX 32
static int _writeback(double *o, const int *tra, const int *trb,
                      const double *tvn, const double *tvl, int nt,
                      const double *drift, int *moved, int nx) {
  int any = 0;
  for (int i = 0; i < nt; i++) if (tvl[i] != tvn[i]) any = 1;
  if (!any) return 0;
  double k0[_WB_MAX];
  int order[_WB_MAX];
  for (int i = 0; i < nt; i++) {
    k0[i] = -fabs(tvl[i] - tvn[i]);
    if (isnan(k0[i])) return -1;
    order[i] = i;
  }
  _sort_by(order, nt, k0, tra, trb);
  int parent[_WB_MAX], innode[_WB_MAX], adjn[_WB_MAX];
  int adjto[_WB_MAX][_WB_MAX];
  double adjs[_WB_MAX][_WB_MAX];
  for (int k = 0; k < nx; k++) { parent[k] = k; innode[k] = 0; adjn[k] = 0; }
  for (int m = 0; m < nt; m++) {
    int i = order[m];
    int ra = tra[i], rb = trb[i];
    double vl = tvl[i];
    int ka = _find(parent, ra), kb = _find(parent, rb);
    if (ka == kb) continue;
    parent[ka] = kb;
    adjto[ra][adjn[ra]] = rb; adjs[ra][adjn[ra]] = -vl; adjn[ra]++;
    adjto[rb][adjn[rb]] = ra; adjs[rb][adjn[rb]] = +vl; adjn[rb]++;
    innode[ra] = 1; innode[rb] = 1;
  }
  int written[_WB_MAX], seen[_WB_MAX], done[_WB_MAX], stack[_WB_MAX];
  for (int k = 0; k < nx; k++) { written[k] = 0; seen[k] = 0; done[k] = 0; }
  /* the components, their members in increasing row order (the set's) */
  for (int r = 0; r < nx; r++) {
    if (!innode[r]) continue;
    int root = _find(parent, r);
    if (done[root]) continue;
    done[root] = 1;
    /* the anchor: min by (0 if pinned else 1, drift, row) */
    int anchor = -1;
    int a0 = 0; double a1 = 0.0;
    for (int n = r; n < nx; n++) {
      if (!innode[n] || _find(parent, n) != root) continue;
      int b0 = moved[n] ? 0 : 1;
      double b1 = drift[n];
      int less;
      if (anchor < 0) less = 1;
      else if (b0 != a0) less = b0 < a0;
      else if (!(b1 == a1)) less = b1 < a1;
      else less = n < anchor;
      if (less) { anchor = n; a0 = b0; a1 = b1; }
    }
    int sp = 0;
    seen[anchor] = 1;
    stack[sp++] = anchor;
    while (sp > 0) {
      int u = stack[--sp];
      for (int e = 0; e < adjn[u]; e++) {
        int v = adjto[u][e];
        double s = adjs[u][e];
        if (seen[v]) continue;
        seen[v] = 1;
        if (moved[v]) continue;
        o[v] = o[u] + s;
        written[v] = 1;
        stack[sp++] = v;
      }
    }
  }
  for (int k = 0; k < nx; k++) if (written[k]) moved[k] = 1;
  return 0;
}
"""


## ---------------------------------------------------------------------------
## The per-class kernel.
## ---------------------------------------------------------------------------

def _c_list(vals):
    return '{' + ', '.join(str(int(v)) for v in vals) + '}'


def render(info):
    """`(csrc, nx, layout)` for the class's limiter, or `Refused`.

    `csrc` is the prelude above and the class's `hdl_fn`; `nx` the length
    of the local state the kernel reads (`i`'s, the class's `len(xsyms)`:
    the rows of the spec index it, and it is longer than the terminal
    count), `layout` the packed parameters' `(n_p, t_index)` -- `i`'s, and
    every parameter function's signature is checked against it."""
    from pycircuit.circuit import hdl
    spec = info.get('limit_spec') or []
    if not spec:
        raise Refused('the class declares no $limit')
    fi = (info.get('funcs') or {}).get('i')
    if fi is None or getattr(fi, '_csrc', None) is None:
        raise Refused('the class carries no C source')
    nx = int(fi._cshape[0])
    n_p, t_index = fi._clayout
    if nx > 32 or len(spec) > 32:
        raise Refused('more rows or probes than the kernel holds')
    stmts, seen = [], set()
    trailing = None
    xnames = None
    pexprs = []                                   # (j, k, expr)
    for j, (rows, kind, move, pfs) in enumerate(spec):
        if kind not in KINDS:
            raise Refused(f'unknown limiter kind {kind!r}')
        if len(pfs) > 2:
            raise Refused('a probe with more than two parameters')
        for k, f in enumerate(pfs):
            par = getattr(f, '_hdl_limit_par', None)
            if par is None:
                raise Refused(f'probe {j} parameter {k} carries no ingredients')
            reach, expr, sig_args, unpack, wants_x = par
            tr = [a.name for a in (sig_args[1:] if wants_x else sig_args)]
            if trailing is None:
                trailing = tr
            elif tr != trailing:
                raise Refused('the parameter functions differ in signature')
            if wants_x:
                xs = [nm for nm, _ in unpack]
                if xnames is None:
                    xnames = xs
                elif xs != xnames:
                    raise Refused('the parameter functions differ in unknowns')
            for sym, e_ in (reach or ()):
                if sym.name not in seen:
                    seen.add(sym.name)
                    stmts.append((sym, e_))
            pexprs.append((j, k, expr))
    if trailing is not None:
        temp = [k for k, nm in enumerate(trailing) if nm == hdl.TEMP.name]
        if len(trailing) != n_p or (temp[0] if temp else None) != t_index:
            raise Refused('the parameter functions are not in the packed layout')
    symmap, array_syms = hdl._c_symmap(stmts, trailing or [], xnames or [])
    printer = hdl._CChainPrinter(symmap, array_syms)
    try:
        lines = [f'  const double L_{sym.name} = {printer.doprint(e_)};'
                 for sym, e_ in stmts]
        pv_lines = [f'  const double pv_{j}_{k} = {printer.doprint(e_)};'
                    for j, k, e_ in pexprs]
    except hdl.CUnsupported as e:
        raise Refused(f'C rendering unavailable: {e}')
    text = '\n'.join(lines + pv_lines)
    if 'tanh(' in text:
        raise Refused('a parameter chain calls tanh (an ulp-level exception)')

    npr = len(spec)
    groups = info.get('limit_groups') or []
    grouped = set()
    for _s, ix in groups:
        grouped.update(ix)
    singles = [i for i in range(npr) if i not in grouped]
    npf = {j: len(pfs) for j, (rows, kind, move, pfs) in enumerate(spec)}
    out = ['/* the hdl limiter of a chained class: x is x0 (the last accepted',
           '   point), out arrives as the copy of the unlimited state and is',
           '   limited in place; 1 declines the call (a NaN among the keys to',
           '   sort: the Python closure answers) */',
           f'int {hdl._C_ENTRY}(const double *x, const double *p, double *out) {{',
           f'  enum {{ NX = {nx}, NP = {npr} }};']
    out.append(text)
    out.append(f'  static const int RA[NP] = {_c_list(r[0][0] for r in spec)};')
    out.append(f'  static const int RB[NP] = {_c_list(r[0][1] for r in spec)};')
    out.append(f'  static const int KIND[NP] = {_c_list(KINDS[r[1]] for r in spec)};')
    out.append(f'  static const int MOVE[NP] = {_c_list(r[2] for r in spec)};')
    p0 = ', '.join(f'pv_{j}_0' if npf[j] > 0 else '0.0' for j in range(npr))
    out.append(f'  const double P0[NP] = {{{p0}}};')
    p1 = ', '.join(f'pv_{j}_1' if npf[j] > 1 else '0.0' for j in range(npr))
    out.append(f'  const double P1[NP] = {{{p1}}};')
    out.append('  double o[NX], drift[NX], k0[NP];')
    out.append('  int moved[NX];')
    out.append('  for (int k = 0; k < NX; k++) {')
    out.append('    o[k] = out[k]; drift[k] = fabs(out[k] - x[k]); moved[k] = 0;')
    out.append('  }')
    for seq, idx in groups:
        ng = len(idx)
        out.append('  {   /* a limit_together group' + (', sequential' if seq else '') + ' */')
        out.append(f'    int order[{ng}] = {_c_list(idx)};')
        if not seq:
            out.append(f'    for (int m = 0; m < {ng}; m++) {{')
            out.append('      int j = order[m];')
            out.append('      double vn = o[RA[j]] - o[RB[j]];')
            out.append('      double vo = x[RA[j]] - x[RB[j]];')
            out.append('      double vl = _lim(KIND[j], vn, vo, P0[j], P1[j]);')
            out.append('      k0[j] = -fabs(vl - vn);')
            out.append('      if (isnan(k0[j])) return 1;')
            out.append('    }')
            out.append(f'    _sort_by(order, {ng}, k0, RA, RB);')
        out.append('    double shift[NX]; int taken[NX];')
        out.append('    for (int k = 0; k < NX; k++) { shift[k] = 0.0; taken[k] = 0; }')
        out.append(f'    int tra[{ng}], trb[{ng}]; double tvn[{ng}], tvl[{ng}]; int nt = 0;')
        out.append(f'    for (int m = 0; m < {ng}; m++) {{')
        out.append('      int j = order[m];')
        out.append('      int ra = RA[j], rb = RB[j];')
        out.append('      double vorig = o[ra] - o[rb];')
        out.append('      double vold = x[ra] - x[rb];')
        out.append('      double vin = (vorig + shift[ra]) - shift[rb];')
        out.append('      double vlim = _lim(KIND[j], vin, vold, P0[j], P1[j]);')
        out.append('      if (vlim != vin) {')
        if seq:
            out.append('        int n = rb;')
        else:
            out.append('        int n = (drift[ra] >= drift[rb]) ? ra : rb;')
            out.append('        if (taken[n]) n = (n == ra) ? rb : ra;')
        out.append('        taken[n] = 1;')
        out.append('        shift[n] = shift[n] + (vlim - vin) * ((n == ra) ? 1.0 : -1.0);')
        out.append('      }')
        out.append('      tra[nt] = ra; trb[nt] = rb; tvn[nt] = vorig; tvl[nt] = vlim; nt++;')
        out.append('    }')
        out.append('    if (_writeback(o, tra, trb, tvn, tvl, nt, drift, moved, NX) < 0) return 1;')
        out.append('  }')
    if singles:
        ns = len(singles)
        out.append('  {   /* the single probes */')
        out.append(f'    int order[{ns}] = {_c_list(singles)};')
        out.append(f'    for (int m = 0; m < {ns}; m++) {{')
        out.append('      int i = order[m];')
        out.append('      double vn = o[RA[i]] - o[RB[i]];')
        out.append('      double vo = x[RA[i]] - x[RB[i]];')
        out.append('      double vl = _lim(KIND[i], vn, vo, P0[i], P1[i]);')
        out.append('      k0[i] = -fabs(vl - vn);')
        out.append('      if (isnan(k0[i])) return 1;')
        out.append('    }')
        out.append(f'    _sort_by(order, {ns}, k0, RA, RB);')
        out.append(f'    for (int m = 0; m < {ns}; m++) {{')
        out.append('      int i = order[m];')
        out.append('      int ra = RA[i], rb = RB[i];')
        out.append('      double vnew = o[ra] - o[rb];')
        out.append('      double vold = x[ra] - x[rb];')
        out.append('      double vlim = _lim(KIND[i], vnew, vold, P0[i], P1[i]);')
        out.append('      if (vlim == vnew) continue;')
        out.append('      int cand = (drift[ra] >= drift[rb]) ? ra : rb;')
        out.append('      if (moved[cand]) cand = (cand == ra) ? rb : ra;')
        out.append('      if (moved[cand]) cand = MOVE[i];')
        out.append('      moved[cand] = 1;')
        out.append('      if (cand == ra) o[ra] = o[rb] + vlim; else o[rb] = o[ra] - vlim;')
        out.append('    }')
        out.append('  }')
    out.append('  for (int k = 0; k < NX; k++) out[k] = o[k];')
    out.append('  return 0;')
    out.append('}')
    return _LIMIT_C + '\n'.join(out) + '\n', nx, (n_p, t_index)


class CLimitKernel(_cb.CKernel):
    """The class's limiter, callable as the closure needs it:
    `kernel(element, x, x0, epar)` -> the limited state as a fresh float64
    array, or None where the closure must answer (a temperature that is
    not one number, a vector that is not one, a declined call)."""

    __slots__ = ()

    def __init__(self, ffi, cfn, nx, layout, key, built_s):
        super().__init__(ffi, cfn, (nx,), nx, layout, key, built_s)

    def __call__(self, element, x, x0, epar):
        d = element.__dict__
        packed = d.get('_hdl_cp')
        if packed is None:
            try:
                packed = self.pack(element)
            except (TypeError, ValueError):
                packed = False
            d['_hdl_cp'] = packed
        if packed is False:
            return None
        p, pcast = packed
        t_index = self.t_index
        if t_index is not None:
            T = getattr(epar, 'T', 300.0)
            ## (a float, an int, a numpy float: one number; the rest --
            ## a 0-d array, an array -- is the closure's)
            if type(T) is not float and type(T) is not int \
                    and type(T) is not np.float64:
                return None
            p[t_index] = T
        if type(x) is np.ndarray and x.dtype is _cb._F64 and x.ndim == 1:
            out = x.copy()
        else:
            try:
                out = np.array(x, dtype=float)
            except (TypeError, ValueError):
                return None
            if out.ndim != 1:
                return None
        if (type(x0) is np.ndarray and x0.dtype is _cb._F64 and x0.ndim == 1
                and x0.flags.c_contiguous):
            x0a = x0
        else:
            try:
                x0a = np.ascontiguousarray(x0, dtype=float)
            except (TypeError, ValueError):
                return None
            if x0a.ndim != 1:
                return None
        nx = self.nx
        if out.shape[0] < nx or x0a.shape[0] < nx:
            return None
        ffi = self.ffi
        if self.cfn(ffi.from_buffer(self.dptr, x0a), pcast,
                    ffi.from_buffer(self.dptr, out)):
            return None
        return out


def source_for(info):
    """A `_build_missing_parallel`-shaped object for the class's limit
    kernel, or None where it is refused (the status is set)."""
    try:
        csrc, nx, layout = render(info)
    except Refused as e:
        info['_c_limit'] = None
        info['_c_limit_status'] = f'numpy ({e})'
        return None
    return types.SimpleNamespace(_csrc=csrc, _nx=nx, _layout=layout)


def bind(cls, info, src=None):
    """Bind the class's limit kernel: `info['_c_limit']` the kernel (or
    None) and `info['_c_limit_status']` 'c' or 'numpy (<why>)'.  A
    limiter that cannot be printed or built refuses ITSELF, never the
    class.  `src` is `source_for(info)`'s object when the caller built
    it already (the parallel build), else it is rendered here."""
    if src is None:
        src = source_for(info)
        if src is None:
            return None
    try:
        ffi, cfn, key, _cold, secs = _cb.load_kernel(src._csrc, limit_cdef())
    except (_cb.CompileError, OSError) as e:
        info['_c_limit'] = None
        info['_c_limit_status'] = f'numpy (compile failed: {e})'
        return None
    kern = CLimitKernel(ffi, cfn, src._nx, src._layout, key, secs)
    info['_c_limit'] = kern
    info['_c_limit_status'] = 'c'
    return kern


def unbind(info):
    """Drop the class's limit kernel: the closure answers again."""
    info['_c_limit'] = None
    info['_c_limit_status'] = 'numpy (detached)'

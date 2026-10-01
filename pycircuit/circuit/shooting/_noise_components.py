"""The noise sources as the noise analyses read them: `CY` sampled on the
orbit or at given states, each leaf element's share of it, the white and
coloured component model, and the per-band classification and roots.

`NoiseComponents(pss, states)` is built per call (`PAC._noise_components`):
the solved `PSS` and the states it samples at (None: the stored orbit) are
its whole input, so nothing is kept on `PAC` between calls.  The helpers it
calls are class attributes as well as module functions, so a test replaces
one on a subclass and hands that subclass to `PAC._noise_components`.

History: `doc/shooting_history.md` (the members' old names on `PAC`).
"""

import numpy as np

from pycircuit.circuit.analysis import remove_row_col
from pycircuit.circuit.circuit import defaultepar
from pycircuit.circuit.simwarnings import CostWarning, ModelWarning, warn

from ._numerics import insert_ref

#: the one key of `NoiseComponents.colour_model`'s components: the whole circuit
JOINT_KEY = ('<circuit>',)


def orbit_states(pss, states=None):
    """Full-width state vectors to sample `CY` at: the stored orbit
    `x(t_k)`, `k = 0..N-1` (the default), or the given `states`."""
    n = pss.cir.n
    irn = pss.irefnode
    if states is None:
        xs = np.asarray(pss.waveform[1], dtype=float)
        nsamp = len(pss.factored_period().steps)
        states = [xs[:, k] for k in range(nsamp)]
    out = []
    for xr in states:
        xr = np.asarray(xr, dtype=float).ravel()
        out.append(xr if xr.shape[0] == n else insert_ref(xr, irn))
    return out


def orbit_midpoints(pss):
    """The full-width state HALFWAY through each step of the stored orbit,
    ``(x(t_k) + x(t_{k+1})) / 2``, `k = 0..N-1` -- where `sign_blind` looks
    besides the samples, so a modulation that changes sign BETWEEN two
    samples is caught when its zero lies near the middle of the step."""
    n = pss.cir.n
    irn = pss.irefnode
    xs = np.asarray(pss.waveform[1], dtype=float)
    nsamp = len(pss.factored_period().steps)
    out = []
    for k in range(nsamp):
        xm = 0.5 * (xs[:, k] + xs[:, (k + 1) % xs.shape[1]])
        out.append(xm if xm.shape[0] == n else insert_ref(xm, irn))
    return out


def colour_fit_frequencies(f, f0):
    ## three to fit, two to verify: one BETWEEN the fit points and one
    ## at the FAR end of the band range the fold reaches (up to
    ## ~(N/2 + lmax) f0), so a shape that is not thermal-plus-flicker
    ## is caught where the model would have been extrapolating
    f = abs(float(f))
    return [2.0 * np.pi * x for x in (max(f, 1e-3 * f0), 3.0 * f0 + f,
                                      10.0 * f0 + f, 2.0 * f0 + f,
                                      150.0 * f0 + f)]


def colour_fit(Cs, ws):
    """`(A, B, EF)` with `C(w) = A + B (w1/w)^EF` entry by entry, from
    `Cs` = samples `(N, n, n)` at the five `ws` of
    `colour_fit_frequencies` (three fit, two verify); None when the
    shape is not thermal-plus-power-law anywhere."""
    from scipy.optimize import brentq
    C1, C2, C3, C4, C5 = Cs
    N, n, _ = C1.shape
    A = np.zeros_like(C1); B = np.zeros_like(C1); EF = np.zeros((N, n, n))
    scale = max(float(np.max(np.abs(C1))), 1e-300)
    w1, w2, w3, w4, w5 = ws
    for k in range(N):
        for i in range(n):
            for j in range(n):
                c1, c2, c3 = C1[k, i, j], C2[k, i, j], C3[k, i, j]
                if abs(c1 - c2) <= 1e-12 * scale and abs(c2 - c3) <= 1e-12 * scale:
                    A[k, i, j] = c1
                    continue
                ratio = (c1 - c2) / (c2 - c3)
                def g(ef, ratio=ratio):
                    g1, g2, g3 = 1.0, (w1 / w2) ** ef, (w1 / w3) ** ef
                    return float(np.real((g1 - g2) / (g2 - g3) - ratio))
                try:
                    ef = brentq(g, 0.05, 4.0, xtol=1e-12)
                except ValueError:
                    return None
                g2 = (w1 / w2) ** ef
                Bv = (c1 - c2) / (1.0 - g2)
                A[k, i, j] = c1 - Bv; B[k, i, j] = Bv; EF[k, i, j] = ef
    for wv, Cv in ((w4, C4), (w5, C5)):
        if float(np.max(np.abs(A + B * (w1 / wv) ** EF - Cv))) > 1e-8 * scale:
            return None
    return A, B, EF


def leaf_cy_stamps(cir, x, w, prefix=(), only=None, epar=defaultepar):
    """Yield `(key, G)`: each LEAF element's `CY(x, w)` stamped into
    `cir`'s full `n x n` space, recursing into sub-circuits.  Their sum
    is `cir.CY(x, w, epar)` (the elements are independent by that method's
    own contract).  `only`: that element's key alone -- nothing else is
    evaluated."""
    n = cir.n
    idx = cir._map_indices_2d
    for inst, el in cir.elements.items():
        if only is not None and tuple(only[:len(prefix) + 1]) != prefix + (inst,):
            continue
        rc = idx.get(inst)
        if rc is None:
            continue
        rows, cols = rc
        subx = np.asarray(x)[cir.elementnodemap[inst]]
        if getattr(el, 'elements', None):
            for key, Gc in leaf_cy_stamps(el, subx, w, prefix + (inst,), only,
                                          epar=epar):
                G = np.zeros((n, n), dtype=complex)
                np.add.at(G, (rows, cols), np.asarray(Gc).ravel())
                yield key, G
        else:
            G = np.zeros((n, n), dtype=complex)
            np.add.at(G, (rows, cols),
                      np.asarray(el.CY(subx, w, epar=epar), dtype=complex).ravel())
            yield prefix + (inst,), G


def leaf_noise_amplitudes(cir, x, w, prefix=(), only=None, epar=defaultepar):
    """Yield `(key, W)`: each leaf element's SIGNED coloured-noise
    amplitudes (`Element.noise_amplitudes`, where it has them) in `cir`'s
    full `n`-row space, `(n, S)`; keyed like `leaf_cy_stamps`."""
    n = cir.n
    for inst, el in cir.elements.items():
        if only is not None and tuple(only[:len(prefix) + 1]) != prefix + (inst,):
            continue
        if cir._map_indices_2d.get(inst) is None:
            continue
        nodemap = np.asarray(cir.elementnodemap[inst])
        subx = np.asarray(x)[nodemap]
        if getattr(el, 'elements', None):
            inner = leaf_noise_amplitudes(el, subx, w, prefix + (inst,), only,
                                          epar=epar)
        else:
            fn = getattr(el, 'noise_amplitudes', None)
            Wc = fn(subx, w, epar=epar) if fn is not None else None
            inner = [] if Wc is None else [(prefix + (inst,), Wc)]
        for key, Wc in inner:
            Wc = np.asarray(Wc, dtype=complex)
            W = np.zeros((n, Wc.shape[1]), dtype=complex)
            np.add.at(W, nodemap, Wc)
            yield key, W


def leaf_access(cir, key, irn):
    """`(element, nodemap, rows, cols, ok)` for a TOP-LEVEL leaf: its
    state slice and where its `CY` lands in the REDUCED matrix (`ok`
    masks the reference row/column out) -- or None (nested, or absent),
    for which `NoiseComponents.one_element_cy` walks the tree."""
    if len(key) != 1 or key[0] not in cir.elements:
        return None
    el = cir.elements[key[0]]
    if getattr(el, 'elements', None):
        return None
    rc = cir._map_indices_2d.get(key[0])
    if rc is None:
        return None
    rows, cols = (np.asarray(a).ravel() for a in rc)
    ok = (rows != irn) & (cols != irn)
    return (el, np.asarray(cir.elementnodemap[key[0]]),
            rows[ok] - (rows[ok] > irn), cols[ok] - (cols[ok] > irn), ok)


def exponent_columns(B, EF, W):
    """`[(W_g, ef_g)]`: the SIGNED columns `W` (at the fit frequency) of a
    component whose power-law exponent differs between its entries,
    grouped by their own exponents -- one column is one fluctuation, so
    one exponent, read at its heaviest entry -- when the groups rebuild
    the component, ``sum_g W_g W_g^H r^ef_g = B r^EF``, at `r = w1/w` of
    10 and 0.1 to 1e-6; else None (a column whose entries mix slopes).
    Each group is then a uniform power law with its own signed columns,
    exactly as a uniform component is taken."""
    W = np.asarray(W, dtype=complex)
    B = np.asarray(B, dtype=complex)
    EF = np.real(np.asarray(EF, dtype=complex))
    groups = []
    for s_ in range(W.shape[2]):
        a2 = np.abs(W[:, :, s_]) ** 2
        if not np.any(a2 > 0.0):
            continue
        k, i = np.unravel_index(int(np.argmax(a2)), a2.shape)
        ef = float(EF[k, i, i])
        for g in groups:
            if abs(g[0] - ef) <= 1e-6 * max(1.0, abs(ef)):
                g[1].append(s_)
                break
        else:
            groups.append((ef, [s_]))
    out = [(W[:, :, idx], ef) for ef, idx in groups]
    for r in (10.0, 0.1):
        lhs = sum(np.einsum('kis,kjs->kij', Wg, Wg.conj()) * r ** ef
                  for Wg, ef in out)
        rhs = B * r ** EF
        if float(np.max(np.abs(lhs - rhs))) > 1e-6 * max(
                float(np.max(np.abs(rhs))), 1e-300):
            return None
    return out


def uniform_exponent(B, EF):
    """The one power-law exponent of a component, or None when its
    non-zero entries carry different ones (then `sqrt(B (w1/w)^EF)` is
    not `(w1/w)^(ef/2) sqrt(B)` and must be taken per band)."""
    ## ⚠ ONLY ENTRIES THAT CARRY WEIGHT VOTE, weighted by what a wrong
    ## exponent COSTS.  The exponent is fitted from differences of `CY`,
    ## so a tiny entry (a MOS flicker source at the sample where Vds
    ## crosses zero) has its exponent in the rounding of the white part
    ## beside it, and that noise goes as 1/weight: no weight cut-off
    ## separates them.  Giving entry i the exponent `ref` misstates the
    ## component by `r_i |(w1/w)^d_i - 1| ~ r_i d_i |ln(w1/w)|` of its
    ## scale (`r` relative weight, `d` deviation); bounded over 50
    ## e-folds of band frequency and held to 1e-9, a genuinely different
    ## exponent (d ~ 1) still fails from a weight of 2e-11 up.  The
    ## reference is the LARGEST entry's exponent, not the first's.
    ## History: `doc/shooting_history.md`, `PAC._uniform_exponent`.
    aB = np.abs(B)
    if not np.any(aB > 0):
        return 0.0
    ref = float(np.real(EF.flat[int(np.argmax(aB))]))
    cost = 50.0 * (aB / float(aB.max())) * np.abs(EF - ref)
    return ref if float(np.max(cost)) <= 1e-9 else None


def split_by_exponent(B, EF):
    """`[(B_g, ef_g)]`: `B` split into the index blocks its nonzero entries
    connect, when each block carries ONE exponent (`uniform_exponent`);
    None when a block mixes exponents (correlated entries of different
    slope, which no split makes independent).  Entries in DISJOINT blocks
    (sources of different slope on branches that share no node) are
    independent components, split exactly.  (`PAC._split_by_exponent` until
    2026-10-01, the covariance family's alone.)"""
    aB = np.max(np.abs(np.asarray(B)), axis=0)
    m = aB.shape[0]
    parent = list(range(m))

    def find(i):
        while parent[i] != i:
            parent[i] = parent[parent[i]]
            i = parent[i]
        return i
    for i, j in zip(*np.nonzero(aB > 0.0)):
        parent[find(i)] = find(j)
    blocks = {}
    for i in range(m):
        if np.any(aB[i] > 0.0) or np.any(aB[:, i] > 0.0):
            blocks.setdefault(find(i), []).append(i)
    out = []
    for idx in blocks.values():
        Bg = np.zeros_like(np.asarray(B))
        ix = np.ix_(range(Bg.shape[0]), idx, idx)
        Bg[ix] = np.asarray(B)[ix]
        EFg = np.zeros_like(np.asarray(EF, dtype=float))
        EFg[ix] = np.asarray(EF, dtype=float)[ix]
        ef = uniform_exponent(Bg, EFg)
        if ef is None:
            return None
        out.append((Bg, float(ef)))
    return out


class ColourBand:
    """One coloured component taken PER BAND FREQUENCY
    (`NoiseComponents.colour_components`): `key`; `kind` -- 'stationary',
    'separable', 'moving' (a per-band element, classified at the probes),
    'density' (a power law whose exponents mix inside one block), 'columns'
    (an element whose white remainder joined the white part: its signed
    columns per band) or 'exact' (unclassified); `root(w)` its columns
    `(K, m, r)` at angular `w`, by the
    shortcut its kind allows; `exact_root(w)` without one; `exact_psd(w)`
    the PSD that root is taken of (None where the columns are the element's
    SIGNED amplitudes); `signed`.  A stationary one carries `supp` (its
    support, for a surface that reads its `CY` directly), a separable one
    `W0` (its columns at f0), `xref`, `pq` and `cref` (its level's
    reference point, entry and value)."""
    __slots__ = ('W0', 'cref', 'exact_psd', 'exact_root', 'key', 'kind',
                 'pq', 'root', 'signed', 'supp', 'xref')

    def __init__(self, key):
        self.key = key
        self.supp = self.W0 = self.xref = self.pq = self.cref = None


class ColourComponents:
    """`NoiseComponents.colour_components`' answer: `fixed` `[(key, W, ef,
    B)]` (a uniform power law through columns `W (K, m, r)`; `B` the PSD
    they are the root of, None where they are signed amplitudes), `bands`
    (`ColourBand`s), `rooted` (the keys factored by a root, sign-tested)
    and `w1` (the power laws' reference frequency)."""
    __slots__ = ('bands', 'fixed', 'rooted', 'w1')

    def __init__(self, w1):
        self.fixed, self.bands, self.rooted, self.w1 = [], [], [], w1


#: `NoiseComponents.perband_classify`'s "derive the mode" default
_DERIVE = object()


def warn_signed_unused(model, where):
    """Warn when an element STATED its signed amplitudes and the fold
    factors that component by sqrt(PSD) anyway: a silent fallback
    reproduces the sign-blind answer, which looks like agreement.  A
    component whose exponent differs between entries is taken by its
    columns grouped per exponent (`exponent_columns`) where they
    rebuild it; only where they do not is it lost.

    History: `doc/shooting_history.md`, `PAC._warn_signed_unused`."""
    signed = getattr(model, 'amplitude', None) or {}
    lost = [key for key, B, EF in (getattr(model, 'flicker', None) or [])
            if key in signed and uniform_exponent(B, EF) is None
            and exponent_columns(B, EF, signed[key]) is None]
    if lost:
        warn(
            '%s: %s states SIGNED coloured-noise amplitudes, but its '
            'power-law exponent is not uniform across its entries and its '
            'columns do not carry one exponent each, so the component is '
            'evaluated per band from sqrt(PSD) -- the SIGN-BLIND fold (the '
            '|m| process).  The result is the pre-2026-09-19 one for this '
            'component, not the signed physics.'
            % (where, ', '.join('.'.join(k) for k in lost)), ModelWarning)


def psd_touches_zero(C):
    """Whether a component's PSD, `C (K, m, m)` over the orbit's points,
    TOUCHES ZERO: a diagonal entry that falls to 1e-2 of its maximum, the
    necessary condition for its modulation to change sign (the pnoise
    fold's test).  Where it does, a factor by the root of the PSD is the
    `|m|` process -- see `warn_sign_blind`."""
    d = np.abs(np.real(np.diagonal(np.asarray(C), axis1=-2, axis2=-1)))
    dmax = d.max(axis=0)
    return bool(np.any((dmax > 0) & (d.min(axis=0) <= 1e-2 * dmax)))


def warn_sign_blind(what, keys):
    """THE one warning for a component factored by the root of its PSD
    whose PSD touches zero (`psd_touches_zero`) and whose element states no
    signed amplitudes: the `|m|` process, exact if the modulation keeps its
    sign, wrong in either direction where it changes sign.  `what` names
    the PAC surface.

    ⚠ Every surface that roots a coloured PSD gives it: the modal and
    lineshape spectra (`colour_groups`), the sample series (from
    2026-09-29; `warn_touch=False` kept it silent) and the covariance
    family (`_coloured_prepare`, from 2026-09-29: its power-law flicker was
    rooted without a word)."""
    warn(
        'PAC.%s: the PSD of %s touches zero along the orbit and the '
        'element states no signed noise amplitudes, so it is factored '
        'by the root of its PSD -- the |m| process: exact if the '
        'modulation keeps its sign, wrong in either direction where it '
        'changes sign.  Only the element knows the sign '
        '(Element.noise_amplitudes).'
        % (what, ', '.join('.'.join(k) for k in keys)), ModelWarning)


def separable(Cs, tol=1e-9):
    """Whether ``C(x_j, w_i) = s_i C(x_j, w_0)`` for every point `j` and
    frequency `i` -- a level that follows the state under a fixed
    spectral shape (`Cs`: one `(K, m, m)` stack per frequency)."""
    C0 = np.asarray(Cs[0], dtype=complex)
    n0 = float(np.vdot(C0, C0).real)
    if n0 == 0.0:
        return False
    for Ci in Cs[1:]:
        Ci = np.asarray(Ci, dtype=complex)
        r = complex(np.vdot(C0, Ci)) / n0
        if float(np.max(np.abs(Ci - r * C0))) > tol * max(
                float(np.max(np.abs(Ci))), 1e-300):
            return False
    return True


def psd_sqrt(Cs):
    """Symmetric PSD square roots of a stack `(..., n, n)`."""
    Cs = np.asarray(Cs, dtype=complex)
    Cs = 0.5 * (Cs + np.conj(np.swapaxes(Cs, -1, -2)))
    lam, U = np.linalg.eigh(Cs)
    return np.einsum('...ik,...k,...jk->...ij', U,
                     np.sqrt(np.clip(np.real(lam), 0.0, None)), U.conj())


def cached_root(cy_at):
    """`w -> psd_sqrt(cy_at(w))`, cached per frequency: the sample
    series reads the same `|f + n f0|` for every instant."""
    cache = {}

    def root(w):
        k = float(w)
        if k not in cache:
            cache[k] = psd_sqrt(cy_at(k))
        return cache[k]
    ## (the caller that knows a frequency is done may clear it:
    ## `_sampled_series`)
    root.cache = cache
    return root


class NoiseComponents(object):
    """The noise sources of `pss`'s circuit at `states` (full-width or
    reduced state vectors; None: the stored orbit samples `x(t_k)`).  See
    the module note."""

    JOINT_KEY = JOINT_KEY
    colour_fit_frequencies = staticmethod(colour_fit_frequencies)
    colour_fit = staticmethod(colour_fit)
    leaf_cy_stamps = staticmethod(leaf_cy_stamps)
    leaf_noise_amplitudes = staticmethod(leaf_noise_amplitudes)
    leaf_access = staticmethod(leaf_access)
    exponent_columns = staticmethod(exponent_columns)
    uniform_exponent = staticmethod(uniform_exponent)
    warn_signed_unused = staticmethod(warn_signed_unused)
    separable = staticmethod(separable)
    split_by_exponent = staticmethod(split_by_exponent)
    psd_sqrt = staticmethod(psd_sqrt)
    cached_root = staticmethod(cached_root)

    def __init__(self, pss, states=None):
        self.pss = pss
        self.states = states
        self.cir = pss.cir
        self.irn = pss.irefnode
        self.keep = np.array([i for i in range(self.cir.n) if i != self.irn])
        self._xs = None

    def at(self, states):
        """The same sources at other `states` (the same class, so a test's
        subclass reaches every sub-evaluation)."""
        return type(self)(self.pss, states)

    @property
    def xs(self):
        """The full-width states sampled (`orbit_states`)."""
        if self._xs is None:
            self._xs = orbit_states(self.pss, self.states)
        return self._xs

    def cy_at_states(self, w):
        """The whole circuit's reduced `CY(x, w)` at the orbit samples or at
        `states`, `(K, m, m)`."""
        Cs = []
        for xf in self.xs:
            cyk = np.asarray(self.cir.CY(xf, w, epar=self.pss.epar), dtype=complex)
            (cyk,) = remove_row_col((cyk,), self.irn, self.pss.toolkit)
            Cs.append(np.asarray(cyk, dtype=complex))
        return np.asarray(Cs, dtype=complex)

    def element_cy_samples(self, w):
        """`{key: (K, m, m)}` -- each NOISY leaf element's reduced `CY(x, w)`
        at the orbit samples `x(t_k)` (the default, indexed like
        `cy_at_states`, which they sum to) or at the given `states`.

        ⚠ AN ELEMENT WHOSE CY IS ZERO AT EVERY SAMPLE IS LEFT OUT (a missing
        key reads as zero): every capacitor, inductor and source stacked a
        dense complex `(K, m, m)` of zeros -- 25 of 26 keys, 146 MB at
        m = 14 (the review's M4, 2026-10-01).  A zero stamp is kept as a
        placeholder until the element shows a non-zero one."""
        keep = self.keep
        out = {}
        for xf in self.xs:
            for key, G in self.leaf_cy_stamps(self.cir, xf, w,
                                              epar=self.pss.epar):
                g = G[np.ix_(keep, keep)]
                out.setdefault(key, []).append(g if np.any(g) else None)
        ## (every leaf, noisy or not, in the circuit's order: `model` keeps
        ## that order, so its sums add in the order they always did)
        self._leaf_order = list(out)
        res = {}
        for key, v in out.items():
            if all(x is None for x in v):
                continue
            shape = next(x for x in v if x is not None).shape
            res[key] = np.asarray([np.zeros(shape, dtype=complex) if x is None
                                   else x for x in v], dtype=complex)
        return res

    def one_element_cy(self, key, w):
        """The reduced `CY(x, w)` of the ONE leaf element `key` at `states`,
        `(K, m, m)` -- nothing else evaluated (a modulated per-band colour
        is read at every point for every band frequency).  A top-level
        leaf is stamped straight into the reduced matrix (`leaf_access`)."""
        acc = self.leaf_access(self.cir, key, self.irn)
        if acc is not None:
            el, nodemap, rr, cc, ok = acc
            m = self.cir.n - 1
            xs = self.xs
            out = np.zeros((len(xs), m, m), dtype=complex)
            for k, xf in enumerate(xs):
                C = np.asarray(el.CY(np.asarray(xf)[nodemap], w,
                                     epar=self.pss.epar),
                               dtype=complex)
                np.add.at(out[k], (rr, cc), C.ravel()[ok])
            return out
        keep = self.keep
        out = []
        for xf in self.xs:
            G = None
            for k, Gk in self.leaf_cy_stamps(self.cir, xf, w, only=key,
                                             epar=self.pss.epar):
                if k == key:
                    G = Gk
            out.append(np.zeros((keep.size, keep.size), dtype=complex) if G is None
                       else G[np.ix_(keep, keep)])
        return np.asarray(out, dtype=complex)

    def one_element_amplitudes(self, key, w):
        """The SIGNED amplitudes of the one element `key` at `states`, `(K, m,
        S)`, or None where it states none."""
        keep = self.keep
        out = []
        for xf in self.xs:
            W = None
            for k, Wk in self.leaf_noise_amplitudes(self.cir, xf, w, only=key,
                                                    epar=self.pss.epar):
                if k == key:
                    W = Wk
            if W is None:
                return None
            out.append(W[keep])
        return np.asarray(out, dtype=complex)

    def signed_amplitudes(self, w1, flicker):
        """`{key: (K, m, S)}` for the flicker components whose element states
        its signed amplitudes AND whose amplitudes rebuild the component:
        `W W^dagger = B` at the fit frequency, to 1e-6.  Anything else keeps
        the square root of its PSD -- the |m| process, warned on as before."""
        keep = self.keep
        acc = {}
        for xf in self.xs:
            for key, W in self.leaf_noise_amplitudes(self.cir, xf, w1,
                                                     epar=self.pss.epar):
                acc.setdefault(key, []).append(W[keep])
        out = {}
        for key, B, _EF in flicker:
            if key not in acc:
                continue
            W = np.asarray(acc[key], dtype=complex)
            if W.shape[0] != B.shape[0] or not np.all(np.isfinite(W)):
                continue
            rebuilt = np.einsum('kis,kjs->kij', W, W.conj())
            if float(np.max(np.abs(rebuilt - B))) <= 1e-6 * float(np.max(np.abs(B))):
                out[key] = W
        return out

    def model(self, f, f0):
        """`CY(x(t), w)` as INDEPENDENT components -- one white and one
        coloured part per leaf element -- for the coloured folds.

        ⚠⚠ Independent sources whose modulations differ do not add under
        a joint square root: `sqrt(A(t) + B)` cross-couples them at
        `(t, t')` (a switch's white noise plus a 1/f source at one node:
        +7.3 % of the total).  One root per component restores
        additivity.

        Returns a callable `w -> (N, m, m)` (the summed `CY`, the contract
        `colour_model` had) carrying `white` (summed), `white_parts`
        `[(key, A)]`, `flicker` `[(key, B, EF)]` (`C = A + B (w1/w)^EF` per
        element, fitted and verified as `colour_fit`), `perband` (keys
        whose colour did not fit: evaluated per band) and `w1`.  ⚠ Within
        one element all white terms share a root, as do all power-law terms
        -- independence is resolved to element x {white, coloured}.  A
        circuit whose `CY` is not the sum of its elements' (an override:
        noise correlated across elements) gets `colour_model`, the same
        interface with the whole circuit as ONE element, or None.

        History: `doc/shooting_history.md`, `PAC._cy_components_model`.
        """
        ws = self.colour_fit_frequencies(f, f0)
        per_w = [self.element_cy_samples(w) for w in ws]
        ## (every element that is noisy at SOME fit frequency, in the order
        ## the circuit lists them; one that is zero at a frequency reads as
        ## zero there)
        keys = [key for key in self._leaf_order
                if any(key in pw for pw in per_w)]
        for pw in per_w:
            ref = next(iter(pw.values()), None)
            for key in keys:
                if key not in pw and ref is not None:
                    pw[key] = np.zeros_like(ref)
        m = self.cir.n - 1
        N = (len(self.pss.factored_period().steps) if self.states is None
             else len(self.states))
        ## ⚠ THE ELEMENTS MUST BE THE CIRCUIT'S CY.  A circuit whose `CY` is
        ## not the sum of its leaf elements' (an override: noise CORRELATED
        ## across elements) cannot be split, and splitting it anyway would
        ## silently analyse a different noise model -- so check at two fit
        ## frequencies and take the whole circuit as one element, saying why.
        for i in (0, len(ws) - 1):
            whole = self.cy_at_states(ws[i])
            parts = sum((per_w[i][key] for key in per_w[i]),
                        np.zeros_like(whole))
            scale = max(float(np.max(np.abs(whole))), 1e-300)
            if float(np.max(np.abs(parts - whole))) > 1e-9 * scale:
                warn(
                    'PAC: this circuit\'s CY is not the sum of its elements\' '
                    '(an override: noise correlated across elements?), so its '
                    'sources cannot be split per element; the whole circuit is '
                    'taken as ONE white and ONE coloured component -- two '
                    'independent COLOURED sources under different modulations '
                    'inside it then do not add (white ones do).', ModelWarning)
                return self.colour_model(f, f0)
        white = np.zeros((N, m, m), dtype=complex)
        white_parts, flicker, perband = [], [], []
        for key in keys:
            fit = self.colour_fit([per_w[i][key] for i in range(len(ws))], ws)
            if fit is None:
                perband.append(key)
                continue
            A, B, EF = fit
            if np.any(A):
                white = white + A
                white_parts.append((key, A))
            if np.any(B):
                flicker.append((key, B, EF))
        ## (an element whose signed amplitudes rebuild its `CY` at every fit
        ## frequency is split by them, one column per fluctuation)
        ## ⚠ EACH PER-BAND ELEMENT'S FACTORING IS DECIDED HERE, ONCE, FOR
        ## EVERY SURFACE (`model.modes`; `perband_mode`).  AND A WHITE
        ## REMAINDER JOINS THE WHITE PART: an element stating the amplitudes
        ## of its coloured sources only (the HDL contract) leaves
        ## ``C - W W^H``, its white sources, which belong with the white
        ## parts over EVERY frequency on every surface -- band-limited, the
        ## covariance read -1.4 % / -1.8 % against a separate white source.
        ## Until 2026-10-01 only the covariance did this; the other surfaces
        ## rooted it per band beside the columns (review O2).  Read at f0,
        ## as the covariance read it.  A remainder that is NOT white (it
        ## moves between the fit frequencies; `perband_mode` warns) stays
        ## beside the columns.
        modes, folded = {}, set()
        w0 = 2.0 * np.pi * float(f0)
        for key in perband:
            Cs_k = [per_w[i][key] for i in range(len(ws))]
            mode = self.perband_mode(key, ws, Cs_k)
            if mode == 'white':
                Rs = []
                for i, w in enumerate(ws):
                    Wi = self.one_element_amplitudes(key, w)
                    Rs.append(np.asarray(Cs_k[i], dtype=complex)
                              - np.einsum('kis,kjs->kij', Wi, Wi.conj()))
                scale = max(float(np.max(np.abs(Rs[0]))), 1e-300)
                if max(float(np.max(np.abs(R_ - Rs[0]))) for R_ in Rs[1:]) \
                        <= 1e-6 * scale:
                    W0 = self.one_element_amplitudes(key, w0)
                    R0 = (self.one_element_cy(key, w0)
                          - np.einsum('kis,kjs->kij', W0, W0.conj()))
                    white = white + R0
                    white_parts.append((key, R0))
                    mode = 'signed'
                    folded.add(key)
            modes[key] = mode
        rooted = [key for key in perband if modes[key] is None]
        if rooted:
            warn(
                'PAC: the noise of %s is not thermal-plus-power-law, so it is '
                'evaluated per band with ONE square root per element: '
                'independent sources INSIDE such an element are not split '
                '(measured 4.2e-4 on an EKV stage, thermal + flicker).'
                % ', '.join('.'.join(k) for k in rooted), CostWarning)
        w1 = ws[0]

        def model(w):
            tot = white.copy()
            ## ⚠ |w|: the stop rule asks at NEGATIVE band frequencies
            ## (f - l f0 < 0), where a non-integer power of a negative base is
            ## NaN and the ratio stop never fires
            ## ⚠ and numpy division: exactly ON a harmonic the folded band is
            ## w = 0, where a Python float division raises ZeroDivisionError
            ## from inside pnoise's harmonic guard instead of letting it see
            ## the non-finite CY and refuse by name
            with np.errstate(divide='ignore', invalid='ignore', over='ignore'):
                ratio = np.float64(w1) / np.abs(np.float64(w))
                for _key, B, EF in flicker:
                    tot = tot + B * ratio ** EF
            if perband:
                ew = self.element_cy_samples(w)
                for key in perband:
                    if key in ew:
                        tot = tot + ew[key]
            return tot
        model.white = white
        model.white_parts = white_parts
        model.flicker = flicker
        model.perband = perband
        model.modes = modes
        model.folded = folded
        ## ⚠ THE SIGN: where the element states its coloured AMPLITUDES, the
        ## folds factor the component with them instead of with sqrt(PSD) --
        ## see `Element.noise_amplitudes` (hdl.py)
        model.amplitude = self.signed_amplitudes(w1, flicker)
        model.w1 = w1
        return model

    def colour_model(self, f, f0):
        """The WHOLE circuit's `CY(x, w)` as ONE white and ONE coloured
        component: `A + B (w1/w)^EF` entry by entry (`colour_fit`, three
        frequencies to fit and two to verify), with the interface of
        `model` under the one key `JOINT_KEY`, at the orbit
        samples or at `states`.  For a circuit whose `CY` is not the sum of
        its elements' -- noise CORRELATED across elements, by an override --
        which cannot be split per element.  None when the whole `CY` is not
        thermal-plus-power-law (the caller then evaluates the circuit per
        band, or refuses).  The exponent is per entry, so a mix of flicker
        exponents across entries is fine.

        ⚠ Independence is resolved to circuit x {white, coloured}.  The white
        part enters linearly wherever it enters (the P-form, Demir's
        functional), so white sources add whatever their modulations; the
        coloured part takes ONE root, so two INDEPENDENT coloured sources
        under different modulations inside such a circuit do not add.

        History: `doc/shooting_history.md`, `PAC._cy_colour_model`."""
        ws = self.colour_fit_frequencies(f, f0)
        fit = self.colour_fit([self.cy_at_states(w) for w in ws], ws)
        if fit is None:
            return None
        A, B, EF = fit
        w1 = ws[0]
        def model(w):
            ## |w| -- see `model`: a negative band frequency
            ## otherwise gives NaN and silently disables pnoise's ratio stop
            with np.errstate(divide='ignore', invalid='ignore', over='ignore'):
                return A + B * (np.float64(w1) / np.abs(np.float64(w))) ** EF
        model.white = A
        model.white_parts = [(self.JOINT_KEY, A)] if np.any(A) else []
        model.flicker = [(self.JOINT_KEY, B, EF)] if np.any(B) else []
        model.perband = []
        model.modes, model.folded = {}, set()
        ## no element states signed amplitudes for a joint root
        model.amplitude = {}
        model.w1 = w1
        return model

    def sign_blind(self, model, keys):
        """The keys among `keys` -- the coloured components a surface
        factors by the ROOT of their PSD -- whose PSD touches zero along the
        orbit (`psd_touches_zero`), read at the orbit's samples AND halfway
        through every step (`orbit_midpoints`), whatever states the surface
        itself samples: one verdict per component and solve, the same on
        every surface (cached on the PSS, `_sign_blind_cache`), warned by
        `warn_sign_blind`.

        What is read is what is rooted: a power-law component's coloured
        part `B` (`colour_fit` at one fixed set of frequencies), a per-band
        element's PSD at f0, the whole circuit's (`JOINT_KEY`; `model` None:
        one root of the whole `CY`) at f0.

        ⚠ A HEURISTIC, NECESSARY NOT SUFFICIENT, AND NO TEST ON THE PSD CAN
        BE BOTH.  The PSD carries `|m|`: `m = t^2` (sign-definite, the root
        exact) and `m = t |t|` (sign-changing, the root O(1) wrong) have ONE
        PSD, so it warns on both.  A switch whose zero falls between the
        points read (narrower than about a quarter step) passes unseen.
        Measured on a Lorentzian into an RC under ``g(V_lo)``, unsigned over
        signed for pnoise(cyclostationary) / sampled / covariance: g = V
        3.44 / 1.050 / 1.040, V|V| 2.98 / 1.033 / 1.025, a tanh switch at
        mid-step 5.06 / 1.117 / 1.101, V^2 and |V| exact on all three --
        wrong or exact alike on every surface, so one verdict serves them
        all (review O2, 2026-10-01; until then three tests, two of them
        missing the V|V| case and a white floor on the same node).

        History: `doc/pss_log_260902.md`, 2026-10-01 (review O2)."""
        keys = list(dict.fromkeys(keys))
        if not keys:
            return []
        pss = self.pss
        cache = getattr(pss, '_sign_blind_cache', None)
        if cache is None:
            cache = {}
            pss._sign_blind_cache = cache
        flick = ({k for k, _B, _EF in model.flicker} if model is not None
                 else set())
        todo = [k for k in keys if (k, k in flick) not in cache]
        if todo:
            f0 = 1.0 / float(pss.period)
            w0 = 2.0 * np.pi * f0
            probe = self.at(orbit_states(pss) + orbit_midpoints(pss))
            ws = self.colour_fit_frequencies(1e-3 * f0, f0)
            for key in todo:
                whole = key == self.JOINT_KEY or model is None
                if key in flick:
                    Cs = [probe.cy_at_states(w) if whole
                          else probe.one_element_cy(key, w) for w in ws]
                    fit = self.colour_fit(Cs, ws)
                    C = fit[1] if fit is not None else Cs[0]
                else:
                    C = (probe.cy_at_states(w0) if whole
                         else probe.one_element_cy(key, w0))
                cache[(key, key in flick)] = psd_touches_zero(C)
        return [k for k in keys if cache[(k, k in flick)]]

    def perband_mode(self, key, ws, Cs=None):
        """How the per-band element `key` is to be factored, decided over
        the frequencies `ws` (`Cs`: its `CY` there, if in hand):

          None      it states no signed amplitudes, or they do not fit its
                    `CY` (warned): the root of its PSD;
          'signed'  its SIGNED amplitudes rebuild its `CY` (`W W^H = C` to
                    1e-6);
          'white'   they rebuild all but a positive semi-definite REMAINDER
                    ``C - W W^H`` -- its sources that state no amplitude, on
                    the HDL contract its WHITE ones (`noise_amplitudes`
                    carries the coloured sources only) -- which is rooted
                    beside the signed columns (`perband_amplitudes`).  A
                    white source has no correlation across the period for a
                    sign to act on, so its root loses nothing.  A remainder
                    that is not white (it moves between `ws`) is sign-blind,
                    and warned.

        ⚠ THE ROOT OF THE PSD IS WRONG TWO WAYS for a modulated per-band
        source: the |m| process where a modulation changes sign, and ONE
        column per point, which merges the element's independent sources.
        Measured on two Lorentzians under two modulations in one element:
        +9 % / +8 % in `sampled_variance` and `covariance` with one
        modulation crossing zero, and still +2.5 % / +5.8 % with both
        positive."""
        modes, Rs = [], []
        for i, w in enumerate(ws):
            W = self.one_element_amplitudes(key, w)
            if W is None:
                return None
            C = (self.one_element_cy(key, w) if Cs is None
                 else np.asarray(Cs[i], dtype=complex))
            R = C - np.einsum('kis,kjs->kij', W, W.conj())
            scale = max(float(np.max(np.abs(C))), 1e-300)
            if float(np.max(np.abs(R))) <= 1e-6 * scale:
                modes.append('signed')
                continue
            Rh = 0.5 * (R + np.conj(np.swapaxes(R, -1, -2)))
            if float(np.min(np.linalg.eigvalsh(Rh))) < -1e-6 * scale:
                warn(
                    'PAC: %s states signed noise amplitudes (Element.'
                    'noise_amplitudes) whose W W^H exceeds its CY, so its PSD '
                    'is rooted instead -- SIGN-BLIND, and its independent '
                    'sources merged.  The two methods disagree.'
                    % '.'.join(key), ModelWarning)
                return None
            modes.append('white')
            Rs.append(R)
        if 'white' not in modes:
            return 'signed'
        if len(Rs) > 1 and max(float(np.max(np.abs(R_ - Rs[0]))) for R_ in Rs[1:]) \
                > 1e-6 * max(float(np.max(np.abs(Rs[0]))), 1e-300):
            warn(
                'PAC: part of the noise of %s states no signed amplitude and '
                'is COLOURED (it changes between band frequencies), so that '
                'part is rooted beside the signed columns -- SIGN-BLIND where '
                'its modulation changes sign (Element.noise_amplitudes '
                'states the sign).' % '.'.join(key), ModelWarning)
        return 'white'

    def perband_amplitudes(self, key, w, mode, C=None):
        """The per-band element `key`'s columns at `w` and `states`, `(K, m,
        r)`, for a `mode` from `perband_mode`: its SIGNED amplitudes (one
        column per independent fluctuation, with the sign of its
        modulation), and for 'white' the root of the remainder ``C - W W^H``
        on the element's support beside them.  `C`: its `CY` at `w`, if in
        hand ('white' only reads it)."""
        W = self.one_element_amplitudes(key, w)
        if mode != 'white':
            return W
        if C is None:
            C = self.one_element_cy(key, w)
        C = np.asarray(C, dtype=complex)
        R = C - np.einsum('kis,kjs->kij', W, W.conj())
        d = np.max(np.abs(np.diagonal(C, axis1=-2, axis2=-1)), axis=0)
        supp = np.nonzero(d > 0.0)[0]
        Z = np.zeros(C.shape[:2] + (supp.size,), dtype=complex)
        Z[:, supp, :] = self.psd_sqrt(R[:, supp][:, :, supp])
        return np.concatenate((np.asarray(W, dtype=complex), Z), axis=-1)

    def perband_classify(self, key, wt, mode=_DERIVE):
        """How the per-band element `key` varies along the orbit, read at
        the probe frequencies `wt`: `(kind, Cs, Ws, mode)`, `Cs` its `CY`
        per probe, `mode` its factoring (`perband_mode`) and `Ws` its
        columns per probe (None for the root of the PSD).

          'stationary'  the same at every point;
          'separable'   the same up to one factor per frequency (a level
                        under a fixed spectral shape);
          'moving'      neither.

        ⚠ FROM THE COLUMNS WHERE STATED: a `CY` that does not move can
        still carry a sign that does (`k(x) = +-1`), and a separable `CY`
        can hold columns whose shapes differ -- either would be read as the
        wrong kind from the `CY`.  `mode`: the factoring already decided
        (`model.modes`), else derived here."""
        Cs = [self.one_element_cy(key, w_) for w_ in wt]
        if mode is _DERIVE:
            mode = self.perband_mode(key, wt, Cs)
        Ws = None if mode is None else [
            self.perband_amplitudes(key, w_, mode, C_)
            for w_, C_ in zip(wt, Cs)]
        X = Cs if Ws is None else Ws
        if all(float(np.max(np.abs(x_ - x_[:1])))
               <= 1e-12 * max(float(np.max(np.abs(x_))), 1e-300) for x_ in X):
            return 'stationary', Cs, Ws, mode
        return ('separable' if self.separable(X) else 'moving'), Cs, Ws, mode

    def perband_root(self, key, mode):
        """`w -> (K, m, r)`: the per-band element `key` at every point for
        every band frequency, cached per frequency -- its columns for a
        `mode` from `perband_mode` (checked at the classification's
        probes; 'signed' is read without the `CY`, which would double the
        cost), else the root of its PSD."""
        cache = {}

        def root(w):
            k = float(w)
            if k not in cache:
                cache[k] = (self.perband_amplitudes(key, k, mode)
                            if mode is not None else self.psd_sqrt(
                                self.one_element_cy(key, k)))
            return cache[k]
        root.cache = cache
        return root

    @staticmethod
    def band_probes(wlo, whi, f0):
        """The (angular) frequencies a per-band element is classified at,
        from a surface's band `[wlo, whi]`: its lower edge, f0, 10 f0, and
        the band's top AND its middle -- a shape that departs only above the
        middle does not read as separable; a probe more can only send a
        source to the exact path.  One set for every surface (review O2,
        2026-10-01; the sample series had no middle, pnoise read f0 alone).
        History: `doc/shooting_history.md`, `PAC._coloured_prepare`."""
        w0 = 2.0 * np.pi * float(f0)
        return sorted({float(wlo), w0, 10.0 * w0, 0.5 * float(whi),
                       float(whi)})

    def _separable_root(self, key, W0, xref, pq, cref):
        """`w -> (K, m, r)`: a separable element's columns at f0 times the
        root of its level at `w` (``C(x, w) = C(x, w_ref) s(w)``), the level
        read at ONE point; cached per frequency (`.cache`)."""
        cache = {}
        pi_, qi = pq
        xr = [xref]

        def root(w):
            k = float(w)
            if k not in cache:
                c = self.at(xr).one_element_cy(key, k)[0, pi_, qi]
                cache[k] = np.sqrt(max(float(np.real(c / cref)), 0.0)) * W0
            return cache[k]
        root.cache = cache
        return root

    def colour_components(self, model, wlo, whi, f0, what, shortcuts=True):
        """THE one grouping of `model`'s COLOURED components, for every
        noise surface (review O2, 2026-10-01: until then three copies -- the
        sample series' / modal / lineshape `colour_groups`, the covariance
        family's `_coloured_prepare` and pnoise's `_cyclostationary_fold` --
        whose probe frequencies, mixed-exponent handling, white remainders
        and sign tests differed).  `[wlo, whi]` is the surface's band
        (angular), which sets the per-band probes (`band_probes`); each
        per-band element's factoring is the one `model` decided
        (`model.modes`).  Returns a `ColourComponents`:

          `fixed`  a uniform power law through FIXED columns: the element's
                   SIGNED amplitudes where it states them (split by their
                   own exponents where they differ, `exponent_columns`),
                   else the root of its PSD (split into the disjoint index
                   blocks of one exponent each, `split_by_exponent`);
          `bands`  the rest, per band frequency (`ColourBand`): a per-band
                   element -- stationary / separable / moving at the probes,
                   or 'exact' with `shortcuts=False` -- or a power law whose
                   exponents mix inside one block ('density', warned for its
                   cost);
          `rooted` the components factored by the ROOT of their PSD, given
                   the one sign-blind verdict (`sign_blind`) and warned here.

        Each surface builds its own numbers from these.  ⚠ An element whose
        white remainder was folded into the white part (`model.folded`) is
        taken through its signed columns alone, per band and without a
        shortcut: its `CY` still holds the remainder, so a level or a
        stationary density read from it would count it twice."""
        signed = getattr(model, 'amplitude', None) or {}
        modes = getattr(model, 'modes', None) or {}
        folded = getattr(model, 'folded', None) or set()
        self.warn_signed_unused(model, 'PAC.%s' % what)
        out = ColourComponents(model.w1)
        rooted = out.rooted
        for key, B, EF in model.flicker:
            ef = self.uniform_exponent(B, EF)
            W = signed.get(key)
            if ef is not None:
                if W is None:
                    rooted.append(key)
                    out.fixed.append((key, np.asarray(self.psd_sqrt(B),
                                                      dtype=complex),
                                      float(ef), B))
                else:
                    out.fixed.append((key, np.asarray(W, dtype=complex),
                                      float(ef), None))
                continue
            split = (self.exponent_columns(B, EF, W) if W is not None
                     else None)
            if split is not None:
                out.fixed.extend((key, np.asarray(Wg, dtype=complex),
                                  float(efg), None) for Wg, efg in split)
                continue
            rooted.append(key)
            parts = self.split_by_exponent(B, EF)
            if parts is not None:
                out.fixed.extend((key, np.asarray(self.psd_sqrt(Bg),
                                                  dtype=complex),
                                  float(efg), Bg) for Bg, efg in parts)
                continue
            band = ColourBand(key)
            band.kind, band.signed = 'density', False
            band.exact_psd = (lambda w, B=B, EF=EF, w1=model.w1:
                              B * (w1 / w) ** EF)
            band.root = band.exact_root = self.cached_root(band.exact_psd)
            out.bands.append(band)
            warn(
                f'PAC.{what}: the coloured noise of {".".join(key)} carries '
                'different power-law exponents in different entries of one '
                'block, so it has no one amplitude to replay; its density is '
                'rooted at every point per band frequency -- costlier.', CostWarning)
        probes = self.band_probes(wlo, whi, f0)
        w0 = 2.0 * np.pi * float(f0)
        sts = self.states if self.states is not None else self.xs
        K = len(sts)
        for key in model.perband:
            mode = modes[key] if key in modes else self.perband_mode(key, probes)
            band = ColourBand(key)
            band.signed = mode is not None
            band.exact_root = self.perband_root(key, mode)
            band.exact_psd = (None if mode is not None else
                              (lambda w, key=key: self.one_element_cy(key, w)))
            if mode is None:
                rooted.append(key)
            if not shortcuts or key in folded:
                band.kind = 'exact' if not shortcuts else 'columns'
                band.root = band.exact_root
                out.bands.append(band)
                continue
            kind, Cs, Ws, _m = self.perband_classify(key, probes, mode=mode)
            band.kind = kind
            i0 = probes.index(w0)
            if kind == 'stationary':
                one = [sts[0]]
                band.root = self.cached_root(
                    lambda w, key=key, one=one, shp=Cs[0].shape[1:]:
                    np.broadcast_to(self.at(one).one_element_cy(key, w),
                                    (K,) + shp))
                C0 = np.asarray(Cs[i0][0], dtype=complex)
                band.supp = np.nonzero(np.any(np.abs(C0) > 0.0, axis=1))[0]
            elif kind == 'separable':
                Cref = Cs[i0]
                band.W0 = self.psd_sqrt(Cref) if Ws is None else Ws[i0]
                jr, pi_, qi = np.unravel_index(int(np.argmax(np.abs(Cref))),
                                               Cref.shape)
                band.xref, band.pq = sts[jr], (pi_, qi)
                band.cref = complex(Cref[jr, pi_, qi])
                band.root = self._separable_root(key, band.W0, band.xref,
                                                 band.pq, band.cref)
            else:
                band.root = band.exact_root
            out.bands.append(band)
        blind = self.sign_blind(model, rooted)
        if blind:
            warn_sign_blind(what, blind)
        return out

    def colour_groups(self, model, wlo, f0, L, what):
        """`colour_components` as the sample series, the modal spectra and
        the lineshape take it: `('fixed', G, s)` -- columns `G (K, m, r)` at
        `states` and a power weight `s(nu)` -- and `('band', root, None)`,
        the columns per band frequency (`ColourBand.root`, cached per
        frequency: the sample series reads the same `|f + n f0|` at every
        instant), over the band up to the grid's `(L + 1/2) f0`.  The FIXED
        groups before the BAND ones is the order the sample series sums
        them in."""
        comps = self.colour_components(
            model, wlo, 2.0 * np.pi * (float(L) + 0.5) * float(f0), f0, what)
        groups = [('fixed', W, lambda nu, ef=ef, w1=comps.w1:
                   (w1 / np.asarray(nu, dtype=float)) ** ef)
                  for _key, W, ef, _B in comps.fixed]
        groups += [('band', band.root, None) for band in comps.bands]
        return groups

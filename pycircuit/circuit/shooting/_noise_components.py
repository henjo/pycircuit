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
import warnings

import numpy as np

from pycircuit.circuit.analysis import remove_row_col

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


def leaf_cy_stamps(cir, x, w, prefix=(), only=None):
    """Yield `(key, G)`: each LEAF element's `CY(x, w)` stamped into
    `cir`'s full `n x n` space, recursing into sub-circuits.  Their sum
    is `cir.CY(x, w)` (the elements are independent by that method's
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
            for key, Gc in leaf_cy_stamps(el, subx, w, prefix + (inst,), only):
                G = np.zeros((n, n), dtype=complex)
                np.add.at(G, (rows, cols), np.asarray(Gc).ravel())
                yield key, G
        else:
            G = np.zeros((n, n), dtype=complex)
            np.add.at(G, (rows, cols),
                      np.asarray(el.CY(subx, w), dtype=complex).ravel())
            yield prefix + (inst,), G


def leaf_noise_amplitudes(cir, x, w, prefix=(), only=None):
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
            inner = leaf_noise_amplitudes(el, subx, w, prefix + (inst,), only)
        else:
            fn = getattr(el, 'noise_amplitudes', None)
            Wc = fn(subx, w) if fn is not None else None
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
        warnings.warn(
            '%s: %s states SIGNED coloured-noise amplitudes, but its '
            'power-law exponent is not uniform across its entries and its '
            'columns do not carry one exponent each, so the component is '
            'evaluated per band from sqrt(PSD) -- the SIGN-BLIND fold (the '
            '|m| process).  The result is the pre-2026-09-19 one for this '
            'component, not the signed physics.'
            % (where, ', '.join('.'.join(k) for k in lost)),
            RuntimeWarning, stacklevel=4)


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
            cyk = np.asarray(self.cir.CY(xf, w), dtype=complex)
            (cyk,) = remove_row_col((cyk,), self.irn, self.pss.toolkit)
            Cs.append(np.asarray(cyk, dtype=complex))
        return np.asarray(Cs, dtype=complex)

    def element_cy_samples(self, w):
        """`{key: (K, m, m)}` -- each leaf element's reduced `CY(x, w)` at the
        orbit samples `x(t_k)` (the default, indexed like `cy_at_states`, which
        they sum to) or at the given `states`."""
        keep = self.keep
        out = {}
        for xf in self.xs:
            for key, G in self.leaf_cy_stamps(self.cir, xf, w):
                out.setdefault(key, []).append(G[np.ix_(keep, keep)])
        return {key: np.asarray(v, dtype=complex) for key, v in out.items()}

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
                C = np.asarray(el.CY(np.asarray(xf)[nodemap], w), dtype=complex)
                np.add.at(out[k], (rr, cc), C.ravel()[ok])
            return out
        keep = self.keep
        out = []
        for xf in self.xs:
            G = None
            for k, Gk in self.leaf_cy_stamps(self.cir, xf, w, only=key):
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
            for k, Wk in self.leaf_noise_amplitudes(self.cir, xf, w, only=key):
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
            for key, W in self.leaf_noise_amplitudes(self.cir, xf, w1):
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
        keys = [key for key in per_w[0]
                if any(np.any(per_w[i][key]) for i in range(len(ws)))]
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
                warnings.warn(
                    'PAC: this circuit\'s CY is not the sum of its elements\' '
                    '(an override: noise correlated across elements?), so its '
                    'sources cannot be split per element; the whole circuit is '
                    'taken as ONE white and ONE coloured component -- two '
                    'independent COLOURED sources under different modulations '
                    'inside it then do not add (white ones do).',
                    RuntimeWarning, stacklevel=3)
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
        rooted = [key for key in perband
                  if self.perband_mode(key, ws,
                                       [per_w[i][key] for i in range(len(ws))])
                  is None]
        if rooted:
            warnings.warn(
                'PAC: the noise of %s is not thermal-plus-power-law, so it is '
                'evaluated per band with ONE square root per element: '
                'independent sources INSIDE such an element are not split '
                '(measured 4.2e-4 on an EKV stage, thermal + flicker).'
                % ', '.join('.'.join(k) for k in rooted),
                RuntimeWarning, stacklevel=3)
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
                    tot = tot + ew[key]
            return tot
        model.white = white
        model.white_parts = white_parts
        model.flicker = flicker
        model.perband = perband
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
        ## no element states signed amplitudes for a joint root
        model.amplitude = {}
        model.w1 = w1
        return model

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
                warnings.warn(
                    'PAC: %s states signed noise amplitudes (Element.'
                    'noise_amplitudes) whose W W^H exceeds its CY, so its PSD '
                    'is rooted instead -- SIGN-BLIND, and its independent '
                    'sources merged.  The two methods disagree.'
                    % '.'.join(key), RuntimeWarning, stacklevel=3)
                return None
            modes.append('white')
            Rs.append(R)
        if 'white' not in modes:
            return 'signed'
        if len(Rs) > 1 and max(float(np.max(np.abs(R_ - Rs[0]))) for R_ in Rs[1:]) \
                > 1e-6 * max(float(np.max(np.abs(Rs[0]))), 1e-300):
            warnings.warn(
                'PAC: part of the noise of %s states no signed amplitude and '
                'is COLOURED (it changes between band frequencies), so that '
                'part is rooted beside the signed columns -- SIGN-BLIND where '
                'its modulation changes sign (Element.noise_amplitudes '
                'states the sign).' % '.'.join(key), RuntimeWarning, stacklevel=3)
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

    def perband_classify(self, key, wt):
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
        wrong kind from the `CY`."""
        Cs = [self.one_element_cy(key, w_) for w_ in wt]
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
        return root

    def perband_root_sampler(self, key, wlo, f0, L):
        """`w -> (K, m, r)`: the COLUMNS of one per-band element at the
        injection points, for the sample series -- the element alone
        (`one_element_cy`), cached per frequency, and classified once
        (`perband_classify`):

          STATIONARY  the same at every point: one point, broadcast;
          SEPARABLE   ``C(x, w) = C(x, w_ref) s(w)``: the columns per point
                      at `w_ref` once, times ``sqrt(s(w))`` from one point;
          otherwise   the element at every point for every frequency.

        The columns are the element's SIGNED amplitudes where it states
        them (`perband_amplitudes`), else the root of its PSD; `.signed`
        on the returned function says which.

        History: `doc/shooting_history.md`, `PAC._perband_root_sampler`."""
        wref = 2.0 * np.pi * f0
        wt = sorted({float(wlo), wref, 20.0 * np.pi * f0,
                     2.0 * np.pi * (float(L) + 0.5) * f0})
        kind, Cs, Ws, mode = self.perband_classify(key, wt)
        states = self.states
        K = len(states)
        if kind == 'stationary':
            one = [states[0]]
            root = self.cached_root(lambda w: np.broadcast_to(
                self.at(one).one_element_cy(key, w), (K,) + Cs[0].shape[1:]))
        elif kind == 'separable':
            Cref = Cs[wt.index(wref)]
            Wref = self.psd_sqrt(Cref) if Ws is None else Ws[wt.index(wref)]
            jr, pi_, qi = np.unravel_index(int(np.argmax(np.abs(Cref))), Cref.shape)
            xr, cref = [states[jr]], complex(Cref[jr, pi_, qi])
            cache = {}

            def root(w):
                k = float(w)
                if k not in cache:
                    c = self.at(xr).one_element_cy(key, k)[0, pi_, qi]
                    cache[k] = np.sqrt(max(float(np.real(c / cref)), 0.0)) * Wref
                return cache[k]
        else:
            root = self.perband_root(key, mode)
        root.signed = Ws is not None
        return root

    def colour_groups(self, model, wlo, f0, L, what, warn_touch=True):
        """The coloured components of `model` as unit processes through
        their own columns: `('fixed', G, s)` -- columns `G (K, m, r)` at
        `states` and a power weight `s(nu)` (a uniform power law, its
        element's SIGNED amplitudes where stated, else the root of its PSD)
        -- or `('band', root, None)`, the columns per band frequency (a
        power law whose exponent varies across its entries; a per-band
        colour, `perband_root_sampler`).  A component factored by the root
        of its PSD whose PSD TOUCHES ZERO along the orbit is warned on: if
        its modulation changes sign there, that root is the `|m|` process
        (`warn_touch`; the sample series has never warned it).  Every band
        root is cached per frequency (`cached_root`): the sample series
        reads the same `|f + n f0|` at every instant.  The FIXED groups
        before the BAND ones is the order the sample series sums them in."""
        signed = getattr(model, 'amplitude', None) or {}
        self.warn_signed_unused(model, 'PAC.%s' % what)

        def touches(C):
            ## the necessary condition for a sign change, as the pnoise fold
            ## asks it: a diagonal entry that falls to 1e-2 of its maximum
            d = np.abs(np.real(np.diagonal(np.asarray(C), axis1=-2, axis2=-1)))
            dmax = d.max(axis=0)
            return bool(np.any((dmax > 0) & (d.min(axis=0) <= 1e-2 * dmax)))
        groups, blind = [], []
        for key, Bc, EF in model.flicker:
            ef = self.uniform_exponent(Bc, EF)
            W = signed.get(key)
            if ef is not None:
                if W is None:
                    if warn_touch and touches(Bc):
                        blind.append(key)
                    W = self.psd_sqrt(Bc)
                groups.append(('fixed', np.asarray(W, dtype=complex),
                               lambda nu, ef=ef, w1=model.w1:
                               (w1 / np.asarray(nu, dtype=float)) ** ef))
            else:
                ## (signed columns grouped by their own exponents, where
                ## they carry one each: a uniform power law per group)
                split = (self.exponent_columns(Bc, EF, W) if W is not None
                         else None)
                if split is not None:
                    for Wg, efg in split:
                        groups.append(('fixed', Wg,
                                       lambda nu, ef=efg, w1=model.w1:
                                       (w1 / np.asarray(nu, dtype=float)) ** ef))
                    continue
                if warn_touch and touches(Bc):
                    blind.append(key)
                groups.append(('band', self.cached_root(
                    lambda w, Bc=Bc, EF=EF, w1=model.w1: Bc * (w1 / w) ** EF),
                    None))
        for key in model.perband:
            root = self.perband_root_sampler(key, wlo, f0, L)
            ## (the test costs a CY evaluation per point: only when it can
            ## warn, and an element's signed amplitudes carry the sign)
            if warn_touch and not root.signed and touches(self.one_element_cy(
                    key, 2.0 * np.pi * f0)):
                blind.append(key)
            groups.append(('band', root, None))
        if blind:
            warnings.warn(
                'PAC.%s: the PSD of %s touches zero along the orbit and the '
                'element states no signed noise amplitudes, so it is factored '
                'by the root of its PSD -- the |m| process: exact if the '
                'modulation keeps its sign, wrong in either direction where it '
                'changes sign.  Only the element knows the sign '
                '(Element.noise_amplitudes).'
                % (what, ', '.join('.'.join(k) for k in blind)),
                RuntimeWarning, stacklevel=4)
        return groups

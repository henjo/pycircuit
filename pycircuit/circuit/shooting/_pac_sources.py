"""The noise sources: CY sampled on the orbit, the per-element and joint colour
models, their fits and roots, the guards, and the period quadrature.
"""
import numpy as np
import warnings
from pycircuit.circuit.analysis import remove_row_col
from ._numerics import periodic_spline_weights


class _NoiseSources(object):
    """The noise sources: CY sampled on the orbit, the per-element and joint
    colour models, their fits and roots, the guards, and the period
    quadrature.  A theme of `PAC` (see `pac.py`)."""

    def _period_dft(self, pss, S):
        """Fourier coefficients over the period of samples `S` `(N, ...)` taken
        at `fp.times[:N]`, laid out like `numpy.fft.fftfreq`.  Uniform grid:
        the index DFT, unchanged.  Non-uniform: the weighted sum at the TRUE
        times with `PSS._period_quadrature`'s weights -- O(N^2), paid only by
        a caller who chose that grid."""
        S = np.asarray(S, dtype=complex)
        fp = pss.factored_period()
        wq = pss._period_quadrature(fp)
        if wq is None or S.shape[0] != len(wq):
            return np.fft.fft(S, axis=0) / S.shape[0]
        N = S.shape[0]
        tms = np.asarray(fp.times, dtype=float)
        ks = np.fft.fftfreq(N, d=1.0 / N)
        E = np.exp(-2j * np.pi * np.outer(ks, (tms[:N] - tms[0]) / float(tms[N] - tms[0]))) * wq[None, :]
        return np.tensordot(E, S, axes=(1, 0))

    #: the one key of `_cy_colour_model`'s components: the whole circuit
    JOINT_KEY = ('<circuit>',)

    def _cy_colour_model(self, pss, f, f0, states=None):
        """The WHOLE circuit's `CY(x, w)` as ONE white and ONE coloured
        component: `A + B (w1/w)^EF` entry by entry (`_colour_fit`, three
        frequencies to fit and two to verify), with the interface of
        `_cy_components_model` under the one key `JOINT_KEY`, at the orbit
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
        ws = self._colour_fit_frequencies(f, f0)
        fit = self._colour_fit([self._cy_at_states(pss, w, states) for w in ws],
                               ws)
        if fit is None:
            return None
        A, B, EF = fit
        w1 = ws[0]
        def model(w):
            ## |w| -- see `_cy_components_model`: a negative band frequency
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

    @staticmethod
    def _colour_fit_frequencies(f, f0):
        ## three to fit, two to verify: one BETWEEN the fit points and one
        ## at the FAR end of the band range the fold reaches (up to
        ## ~(N/2 + lmax) f0), so a shape that is not thermal-plus-flicker
        ## is caught where the model would have been extrapolating
        f = abs(float(f))
        return [2.0 * np.pi * x for x in (max(f, 1e-3 * f0), 3.0 * f0 + f,
                                          10.0 * f0 + f, 2.0 * f0 + f,
                                          150.0 * f0 + f)]

    @staticmethod
    def _colour_fit(Cs, ws):
        """`(A, B, EF)` with `C(w) = A + B (w1/w)^EF` entry by entry, from
        `Cs` = samples `(N, n, n)` at the five `ws` of
        `_colour_fit_frequencies` (three fit, two verify); None when the
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

    @classmethod
    def _leaf_cy_stamps(cls, cir, x, w, prefix=(), only=None):
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
                for key, Gc in cls._leaf_cy_stamps(el, subx, w, prefix + (inst,),
                                                   only):
                    G = np.zeros((n, n), dtype=complex)
                    np.add.at(G, (rows, cols), np.asarray(Gc).ravel())
                    yield key, G
            else:
                G = np.zeros((n, n), dtype=complex)
                np.add.at(G, (rows, cols),
                          np.asarray(el.CY(subx, w), dtype=complex).ravel())
                yield prefix + (inst,), G

    @classmethod
    def _leaf_noise_amplitudes(cls, cir, x, w, prefix=(), only=None):
        """Yield `(key, W)`: each leaf element's SIGNED coloured-noise
        amplitudes (`Element.noise_amplitudes`, where it has them) in `cir`'s
        full `n`-row space, `(n, S)`; keyed like `_leaf_cy_stamps`."""
        n = cir.n
        for inst, el in cir.elements.items():
            if only is not None and tuple(only[:len(prefix) + 1]) != prefix + (inst,):
                continue
            if cir._map_indices_2d.get(inst) is None:
                continue
            nodemap = np.asarray(cir.elementnodemap[inst])
            subx = np.asarray(x)[nodemap]
            if getattr(el, 'elements', None):
                inner = cls._leaf_noise_amplitudes(el, subx, w, prefix + (inst,),
                                                   only)
            else:
                fn = getattr(el, 'noise_amplitudes', None)
                Wc = fn(subx, w) if fn is not None else None
                inner = [] if Wc is None else [(prefix + (inst,), Wc)]
            for key, Wc in inner:
                Wc = np.asarray(Wc, dtype=complex)
                W = np.zeros((n, Wc.shape[1]), dtype=complex)
                np.add.at(W, nodemap, Wc)
                yield key, W

    def _signed_amplitudes(self, pss, w1, flicker, states=None):
        """`{key: (K, m, S)}` for the flicker components whose element states
        its signed amplitudes AND whose amplitudes rebuild the component:
        `W W^dagger = B` at the fit frequency, to 1e-6.  Anything else keeps
        the square root of its PSD -- the |m| process, warned on as before."""
        irn = pss.irefnode
        keep = np.array([i for i in range(pss.cir.n) if i != irn])
        acc = {}
        for xf in self._orbit_states(pss, states):
            for key, W in self._leaf_noise_amplitudes(pss.cir, xf, w1):
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

    @staticmethod
    def _exponent_columns(B, EF, W):
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

    @classmethod
    def _warn_signed_unused(cls, model, where):
        """Warn when an element STATED its signed amplitudes and the fold
        factors that component by sqrt(PSD) anyway: a silent fallback
        reproduces the sign-blind answer, which looks like agreement.  A
        component whose exponent differs between entries is taken by its
        columns grouped per exponent (`_exponent_columns`) where they
        rebuild it; only where they do not is it lost.

        History: `doc/shooting_history.md`, `PAC._warn_signed_unused`."""
        signed = getattr(model, 'amplitude', None) or {}
        lost = [key for key, B, EF in (getattr(model, 'flicker', None) or [])
                if key in signed and cls._uniform_exponent(B, EF) is None
                and cls._exponent_columns(B, EF, signed[key]) is None]
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

    @staticmethod
    def _orbit_states(pss, states=None):
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
            out.append(xr if xr.shape[0] == n
                       else np.concatenate((xr[:irn], np.zeros(1), xr[irn:])))
        return out

    @staticmethod
    def _leaf_access(cir, key, irn):
        """`(element, nodemap, rows, cols, ok)` for a TOP-LEVEL leaf: its
        state slice and where its `CY` lands in the REDUCED matrix (`ok`
        masks the reference row/column out) -- or None (nested, or absent),
        for which `_one_element_cy` walks the tree."""
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

    def _one_element_cy(self, pss, key, w, states):
        """The reduced `CY(x, w)` of the ONE leaf element `key` at `states`,
        `(K, m, m)` -- nothing else evaluated (a modulated per-band colour
        is read at every point for every band frequency).  A top-level
        leaf is stamped straight into the reduced matrix (`_leaf_access`)."""
        irn = pss.irefnode
        acc = self._leaf_access(pss.cir, key, irn)
        if acc is not None:
            el, nodemap, rr, cc, ok = acc
            m = pss.cir.n - 1
            xs = self._orbit_states(pss, states)
            out = np.zeros((len(xs), m, m), dtype=complex)
            for k, xf in enumerate(xs):
                C = np.asarray(el.CY(np.asarray(xf)[nodemap], w), dtype=complex)
                np.add.at(out[k], (rr, cc), C.ravel()[ok])
            return out
        keep = np.array([i for i in range(pss.cir.n) if i != irn])
        out = []
        for xf in self._orbit_states(pss, states):
            G = None
            for k, Gk in self._leaf_cy_stamps(pss.cir, xf, w, only=key):
                if k == key:
                    G = Gk
            out.append(np.zeros((keep.size, keep.size), dtype=complex) if G is None
                       else G[np.ix_(keep, keep)])
        return np.asarray(out, dtype=complex)

    def _one_element_amplitudes(self, pss, key, w, states):
        """The SIGNED amplitudes of the one element `key` at `states`, `(K, m,
        S)`, or None where it states none."""
        irn = pss.irefnode
        keep = np.array([i for i in range(pss.cir.n) if i != irn])
        out = []
        for xf in self._orbit_states(pss, states):
            W = None
            for k, Wk in self._leaf_noise_amplitudes(pss.cir, xf, w, only=key):
                if k == key:
                    W = Wk
            if W is None:
                return None
            out.append(W[keep])
        return np.asarray(out, dtype=complex)

    def _perband_mode(self, pss, key, ws, states, Cs=None):
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
                    beside the signed columns (`_perband_amplitudes`).  A
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
            W = self._one_element_amplitudes(pss, key, w, states)
            if W is None:
                return None
            C = (self._one_element_cy(pss, key, w, states) if Cs is None
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

    def _perband_amplitudes(self, pss, key, w, states, mode, C=None):
        """The per-band element `key`'s columns at `w` and `states`, `(K, m,
        r)`, for a `mode` from `_perband_mode`: its SIGNED amplitudes (one
        column per independent fluctuation, with the sign of its
        modulation), and for 'white' the root of the remainder ``C - W W^H``
        on the element's support beside them.  `C`: its `CY` at `w`, if in
        hand ('white' only reads it)."""
        W = self._one_element_amplitudes(pss, key, w, states)
        if mode != 'white':
            return W
        if C is None:
            C = self._one_element_cy(pss, key, w, states)
        C = np.asarray(C, dtype=complex)
        R = C - np.einsum('kis,kjs->kij', W, W.conj())
        d = np.max(np.abs(np.diagonal(C, axis1=-2, axis2=-1)), axis=0)
        supp = np.nonzero(d > 0.0)[0]
        Z = np.zeros(C.shape[:2] + (supp.size,), dtype=complex)
        Z[:, supp, :] = self._psd_sqrt(R[:, supp][:, :, supp])
        return np.concatenate((np.asarray(W, dtype=complex), Z), axis=-1)

    def _perband_classify(self, pss, key, states, wt):
        """How the per-band element `key` varies along the orbit, read at
        the probe frequencies `wt`: `(kind, Cs, Ws, mode)`, `Cs` its `CY`
        per probe, `mode` its factoring (`_perband_mode`) and `Ws` its
        columns per probe (None for the root of the PSD).

          'stationary'  the same at every point;
          'separable'   the same up to one factor per frequency (a level
                        under a fixed spectral shape);
          'moving'      neither.

        ⚠ FROM THE COLUMNS WHERE STATED: a `CY` that does not move can
        still carry a sign that does (`k(x) = +-1`), and a separable `CY`
        can hold columns whose shapes differ -- either would be read as the
        wrong kind from the `CY`."""
        Cs = [self._one_element_cy(pss, key, w_, states) for w_ in wt]
        mode = self._perband_mode(pss, key, wt, states, Cs)
        Ws = None if mode is None else [
            self._perband_amplitudes(pss, key, w_, states, mode, C_)
            for w_, C_ in zip(wt, Cs)]
        X = Cs if Ws is None else Ws
        if all(float(np.max(np.abs(x_ - x_[:1])))
               <= 1e-12 * max(float(np.max(np.abs(x_))), 1e-300) for x_ in X):
            return 'stationary', Cs, Ws, mode
        return ('separable' if self._separable(X) else 'moving'), Cs, Ws, mode

    def _perband_root(self, pss, key, states, mode):
        """`w -> (K, m, r)`: the per-band element `key` at every point for
        every band frequency, cached per frequency -- its columns for a
        `mode` from `_perband_mode` (checked at the classification's
        probes; 'signed' is read without the `CY`, which would double the
        cost), else the root of its PSD."""
        cache = {}

        def root(w):
            k = float(w)
            if k not in cache:
                cache[k] = (self._perband_amplitudes(pss, key, k, states, mode)
                            if mode is not None else self._psd_sqrt(
                                self._one_element_cy(pss, key, k, states)))
            return cache[k]
        return root

    @staticmethod
    def _separable(Cs, tol=1e-9):
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

    def _element_cy_samples(self, pss, w, states=None):
        """`{key: (K, m, m)}` -- each leaf element's reduced `CY(x, w)` at the
        orbit samples `x(t_k)` (the default, indexed like `_cy_at_states`, which
        they sum to) or at the given `states`."""
        irn = pss.irefnode
        n = pss.cir.n
        keep = np.array([i for i in range(n) if i != irn])
        out = {}
        for xf in self._orbit_states(pss, states):
            for key, G in self._leaf_cy_stamps(pss.cir, xf, w):
                out.setdefault(key, []).append(G[np.ix_(keep, keep)])
        return {key: np.asarray(v, dtype=complex) for key, v in out.items()}

    def _cy_at_states(self, pss, w, states=None):
        """The whole circuit's reduced `CY(x, w)` at the orbit samples or at
        `states`, `(K, m, m)`."""
        irn = pss.irefnode
        Cs = []
        for xf in self._orbit_states(pss, states):
            cyk = np.asarray(pss.cir.CY(xf, w), dtype=complex)
            (cyk,) = remove_row_col((cyk,), irn, pss.toolkit)
            Cs.append(np.asarray(cyk, dtype=complex))
        return np.asarray(Cs, dtype=complex)

    def _cy_components_model(self, pss, f, f0, states=None):
        """`CY(x(t), w)` as INDEPENDENT components -- one white and one
        coloured part per leaf element -- for the coloured folds.

        ⚠⚠ Independent sources whose modulations differ do not add under
        a joint square root: `sqrt(A(t) + B)` cross-couples them at
        `(t, t')` (a switch's white noise plus a 1/f source at one node:
        +7.3 % of the total).  One root per component restores
        additivity.

        Returns a callable `w -> (N, m, m)` (the summed `CY`, the contract
        `_cy_colour_model` had) carrying `white` (summed), `white_parts`
        `[(key, A)]`, `flicker` `[(key, B, EF)]` (`C = A + B (w1/w)^EF` per
        element, fitted and verified as `_colour_fit`), `perband` (keys
        whose colour did not fit: evaluated per band) and `w1`.  ⚠ Within
        one element all white terms share a root, as do all power-law terms
        -- independence is resolved to element x {white, coloured}.  A
        circuit whose `CY` is not the sum of its elements' (an override:
        noise correlated across elements) gets `_cy_colour_model`, the same
        interface with the whole circuit as ONE element, or None.

        History: `doc/shooting_history.md`, `PAC._cy_components_model`.
        """
        ws = self._colour_fit_frequencies(f, f0)
        per_w = [self._element_cy_samples(pss, w, states) for w in ws]
        keys = [key for key in per_w[0]
                if any(np.any(per_w[i][key]) for i in range(len(ws)))]
        m = pss.cir.n - 1
        N = (len(pss.factored_period().steps) if states is None
             else len(states))
        ## ⚠ THE ELEMENTS MUST BE THE CIRCUIT'S CY.  A circuit whose `CY` is
        ## not the sum of its leaf elements' (an override: noise CORRELATED
        ## across elements) cannot be split, and splitting it anyway would
        ## silently analyse a different noise model -- so check at two fit
        ## frequencies and take the whole circuit as one element, saying why.
        for i in (0, len(ws) - 1):
            whole = self._cy_at_states(pss, ws[i], states)
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
                return self._cy_colour_model(pss, f, f0, states)
        white = np.zeros((N, m, m), dtype=complex)
        white_parts, flicker, perband = [], [], []
        for key in keys:
            fit = self._colour_fit([per_w[i][key] for i in range(len(ws))], ws)
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
                  if self._perband_mode(pss, key, ws, states,
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
                ew = self._element_cy_samples(pss, w, states)
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
        model.amplitude = self._signed_amplitudes(pss, w1, flicker, states)
        model.w1 = w1
        return model

    @staticmethod
    def _uniform_exponent(B, EF):
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

    def _cy_cycle_averaged(self, pss, w):
        """`CY` time-averaged over the orbit — Hull & Meyer's construction.

        ⚠⚠ VALID FOR GENTLE MODULATION ONLY, AND IT FAILS AS A FACTOR, NOT A
        PERCENTAGE: on a series switch + shunt capacitor it reads 16x a
        reference simulator at `goff/gon` = 1e-6 (1.000 unmodulated).  The
        averaged source injects `4kT <g>` for the WHOLE period, including
        the hold phase, where Hull & Meyer's own condition fails by six
        orders.  So this route is for a mixer's `gm` or a bias-dependent
        shot noise -- not for a switch.  A switch's noise is reachable
        exactly: `covariance` and `oscillator_covariance` evaluate `CY` at
        every step and need no averaging.

        This is what `_cy_reduced` refuses, done instead of refused, after
        Hull & Meyer (1993): *"cyclostationary noise sources, such as shot
        noise, may be modeled as MODULATED STATIONARY NOISE SOURCES.  The
        impulse response that is calculated INCLUDES THE EFFECT OF THIS
        MODULATION"*, with the stationary source at the cycle-averaged
        value.  The modulation is carried by the RESPONSE, `H_l`, rather
        than by the sources: ONE source per device.

        ⚠ Its condition: *"NONE OF THE LARGE-SIGNAL STATE VARIABLES MAY
        CHANGE SIGNIFICANTLY OVER THE DECAY TIME OF THE IMPULSE
        RESPONSE."*  A ringing (high-Q) impulse response breaks it, so it
        degrades as `lambda_2 -> 1`; `info` reports `|lambda_2|` so the
        caller can see which regime they are in.

        ⚠ SAMPLED ON THE ORBIT, NOT AT THE OPERATING POINT: `CY` is
        evaluated at every stored state and averaged with the step weights
        -- the same quadrature `diffusion_constant` uses, so the two remain
        comparable.

        History: `doc/shooting_history.md`, `PAC._cy_cycle_averaged`.
        """
        irn = pss.irefnode
        fp = pss.factored_period()
        tms = np.asarray(fp.times, dtype=float)
        T = float(fp.T)
        xs = np.asarray(pss.waveform[1], dtype=float)
        nsamp = min(len(tms) - 1, xs.shape[1])
        hs = self._period_weights(tms, nsamp, T, pss)
        acc = None
        for k in range(nsamp):
            ## ⚠ `waveform` is FULL width: `xs[:m]` then a second zero reads
            ## every unknown past the reference one slot late
            ## History: `doc/shooting_history.md`, `PAC._cy_cycle_averaged`.
            xf = self._orbit_states(pss, [xs[:, k]])[0]
            cyk = np.asarray(pss.cir.CY(xf, w), dtype=complex)
            (cyk,) = remove_row_col((cyk,), irn, pss.toolkit)
            cyk = np.asarray(cyk, dtype=complex) * hs[k]
            acc = cyk if acc is None else acc + cyk
        if acc is None:
            raise NotImplementedError(
                'PAC: the orbit has no stored samples to average CY over.')
        return acc / float(hs[:nsamp].sum())

    def _cy_reduced(self, pss, w):
        """`CY` with the reference node removed, refusing a moving one.

        ⚠ THE CHECK IS THE POINT.  A bias-dependent `CY` makes the sources
        CYCLOSTATIONARY, and then sidebands stop adding in power -- the
        stationary sum this class computes would be the wrong model, not
        merely an inaccurate one, and nothing downstream would say so.
        Sampled at three states on the converged orbit rather than argued
        from the element types, because a compact model's `CY` reads `x`
        and the discrete library's does not.

        ⚠⚠ IN EFFECT THIS REFUSES MOS pnoise: no physically correct MOS
        noise model has a state-independent `CY` (thermal `4kT gamma g_d0`
        with `g_d0` bias-dependent, flicker `I_D^AF`, gate shot `2qI_G`,
        trap rates reading the terminal voltages -- Mahmutoglu & Demir
        2015).  `PspMosLongChannel` is noiseless only by default (`fnt =
        0`).  The routes past the refusal are `cyclostationary=True` (exact
        for thermal and shot noise; flicker is the `|m|` fold, see
        `pnoise`) and `modulated=True` (Hull & Meyer's cycle average, see
        `_cy_cycle_averaged`).

        ⚠ This same check keeps the Ito/Stratonovich choice out of reach
        (`CY = GG^T`, so a state-dependent `CY` is a state-dependent `G`).
        Past it the two interpretations diverge, and Demir's escape -- "the
        noise signals are small compared with the deterministic signals" --
        may NOT carry for trap noise (a two-state Markov chain, not a small
        perturbation).  The tell would be a discrepancy in a MEAN but not in
        a variance.  (The surfaces that read `CY` state by state --
        `covariance`, `pnoise(cyclostationary=True)`, `diffusion_constant`,
        `modal_spectrum` -- report DIFFUSION, where the two agree, and no
        mean.)

        History: `doc/shooting_history.md`, `PAC._cy_reduced`.
        """
        irn = pss.irefnode
        fp = pss.factored_period()
        ## ⚠ THREE STATES ON THE ORBIT: the stored state half a period in,
        ## `x_last` and `x_prev`.  Never the zero vector -- it is on the
        ## orbit only by accident, and a switch model reading `goff` at
        ## v(ck) = 0 would refuse a linear time-invariant RC held by a DC
        ## clock as cyclostationary.
        _W = np.asarray(pss.waveform[1], dtype=float)
        _mid = np.delete(_W[:, _W.shape[1] // 2], irn, axis=0)
        states = [np.asarray(_mid, dtype=float).ravel()[:pss.cir.n - 1],
                  np.asarray(fp.x_last, dtype=float).ravel(),
                  np.asarray(fp.x_prev, dtype=float).ravel()[:pss.cir.n - 1]]
        mats = []
        for xr in states:
            xf = np.concatenate((xr[:irn], np.zeros(1), xr[irn:]))
            cy = np.asarray(pss.cir.CY(xf, w), dtype=complex)
            (cy,) = remove_row_col((cy,), irn, pss.toolkit)
            mats.append(np.asarray(cy, dtype=complex))
        scale = max(float(np.max(np.abs(mats[0]))), 1e-300)
        for other in mats[1:]:
            drift = float(np.max(np.abs(other - mats[0]))) / scale
            if drift > 1e-9:
                raise NotImplementedError(
                    'PAC: this circuit has a BIAS-DEPENDENT CY (varies by '
                    '%.3g over the orbit), so its noise sources are '
                    'cyclostationary. The sidebands are then correlated '
                    'through the window Fourier coefficients and no longer '
                    'add in power, so the stationary sum here would be the '
                    'wrong model rather than an imprecise one. '
                    '⚠ THIS IS THE NORMAL CASE FOR ANY COMPACT MOS MODEL, '
                    'not an exotic one: thermal channel noise is '
                    '4kT.gamma.gd0 with gd0 bias-dependent, flicker goes '
                    'as I_D^AF, gate shot noise as 2qI_G, and trap capture '
                    'and emission rates read the terminal voltages. So '
                    'this is in effect a refusal of MOS pnoise, and the '
                    'route out is the CYCLOSTATIONARY construction rather '
                    'than a different device model -- BUILT: pass '
                    'cyclostationary=True (the PSD\'s own harmonics, exact; '
                    'modulated=True is the cycle average, which drops the '
                    'sideband correlation and read 0.53 of the truth on a '
                    'driven multiplier). Hull & Meyer (1993) '
                    'make it affordable -- ONE stationary source per '
                    'device at the CYCLE-AVERAGED current, with the '
                    'modulation carried by the impulse response, valid '
                    'while no large-signal state variable changes much '
                    'over the impulse response decay time. ⚠ THAT ROUTE '
                    'IS BUILT: pass modulated=True to use it. It is a '
                    'MODEL CHOICE with the validity condition above, not '
                    'a tolerance relaxation, which is why it is opt-in '
                    'and why this refusal is the default.'
                    % drift)
        return mats[0]

    ## How close to a harmonic of `f0` counts as "on" it, as a fraction of
    ## `f0` (`HARMONIC_GUARD`, below).  The deflated solve's conditioning is
    ## FLAT down to 1e-9 of `f0`, so the guard excludes only what has no
    ## finite answer: at an EXACT harmonic `1/(1 - alpha)` is a division by
    ## zero and the response is genuinely unbounded.
    ## History: `doc/shooting_history.md`, `PAC.HARMONIC_GUARD`.
    def _cy_at(self, pss, w, xr):
        """`CY` at ONE reduced state `xr`, reference row/column removed, with
        NO cyclostationarity check -- for the routes that evaluate the
        source at every step and therefore model a modulated source
        exactly (`covariance`, `oscillator_covariance`).  The stationary
        sum in `pnoise` cannot, which is why `_cy_reduced` refuses there.
        """
        irn = pss.irefnode
        xr = np.asarray(xr, dtype=float).ravel()[:pss.cir.n - 1]
        xf = np.concatenate((xr[:irn], np.zeros(1), xr[irn:]))
        cy = np.asarray(pss.cir.CY(xf, w), dtype=complex)
        (cy,) = remove_row_col((cy,), irn, pss.toolkit)
        return np.asarray(cy, dtype=complex)

    def _refuse_coloured(self, pss, what, instead=None):
        """Refuse a coloured source where the machinery assumes WHITE.
        (`covariance` and `event_jitter` take a band instead and do not
        come here with one -- `_coloured_prepare`.)

        ⚠ THE TRAP IS THAT NOTHING ELSE WOULD OBJECT. The Lyapunov
        recursion, `diffusion_constant` and eq (22)'s collapse all read
        `CY` at ONE frequency and treat it as the noise intensity at every
        frequency; a coloured source folded that way returns a plausible
        number, not an error (A4d names exactly this shape). Detected by
        evaluating the reduced `CY` at two frequencies -- colour is
        frequency dependence, bias dependence is what `_cy_reduced`
        refuses separately.
        """
        if self._coloured_present(pss):
            raise NotImplementedError(
                'PAC.%s: a noise source in this circuit is COLOURED (its CY '
                'differs between w0 and 10 w0), and this routine assumes '
                'white sources -- it would fold CY at one frequency as if '
                'it held at every frequency and return a plausible wrong '
                'number. Use the frequency-resolved surfaces (pnoise, '
                'sampled_variance, phase_psd/coloured_diffusion on an '
                'oscillator, where 1/f noise makes the phase growth '
                'non-diffusive), covariance/event_jitter with a band '
                '(fmin, fmax) on a driven circuit, or the white-through-filter '
                'form of the source.%s' % (what, (' Here: ' + instead)
                                            if instead else ''))

    def _coloured_present(self, pss):
        """Whether any noise source of the circuit is COLOURED -- see
        `_refuse_coloured`, which asks this."""
        w1 = 2.0 * np.pi / float(pss.period)
        ## ⚠ Colour is asked at fixed state, two frequencies: it is separable
        ## from the bias question, and asking it through `_cy_reduced` would
        ## refuse every MODULATED source the covariance routes (which
        ## evaluate `CY` per step) handle exactly.
        ## ⚠⚠ PER ENTRY, AND AT MORE THAN ONE STATE: each entry is judged
        ## against ITS OWN magnitude (a drain's flicker colour sits orders
        ## below a gate resistor's white 4kT/rg in the same matrix), at
        ## `x_last` and at states spread over the orbit (a switch OFF at
        ## t = 0 is white there and coloured elsewhere).  Exact zeros and
        ## white entries give identical values at both frequencies, so
        ## neither can fire.
        ## `_cy_at` takes a REDUCED state: the reference row is removed here.
        ## History: `doc/shooting_history.md`, `PAC._refuse_coloured`.
        _xl = np.asarray(pss.factored_period().x_last, dtype=float).ravel()
        _states = [_xl]
        _wf = getattr(pss, 'waveform', None)
        if _wf is not None:
            _W = np.delete(np.asarray(_wf[1], dtype=float), pss.irefnode,
                           axis=0)
            for _k in sorted(set(np.linspace(0, _W.shape[1] - 1,
                                             8).astype(int))):
                _states.append(_W[:, _k])
        for _xr in _states:
            c1 = self._cy_at(pss, w1, _xr)
            c2 = self._cy_at(pss, 10.0 * w1, _xr)
            den = np.maximum(np.abs(c1), np.abs(c2))
            if np.any(np.abs(c1 - c2) > 1e-9 * den):
                return True
        return False

    @classmethod
    def _split_by_exponent(cls, B, EF):
        """`[(B_g, ef_g)]`: `B` split into the index blocks its nonzero entries
        connect, when each block carries ONE exponent (`_uniform_exponent`);
        None when a block mixes exponents (correlated entries of different
        slope, which no split makes independent)."""
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
            ef = cls._uniform_exponent(Bg, EFg)
            if ef is None:
                return None
            out.append((Bg, float(ef)))
        return out

    @classmethod
    def _power_law_weights(cls, nus, ef, richardson=True):
        """Weights `q_i` with ``int nu^-ef g(nu) dnu = sum_i q_i g(nu_i)``
        EXACT for `g` linear in ``ln nu`` between the points -- the power law
        integrated analytically, the response interpolated.  On ``[u_i,
        u_i + h]`` in ``u = ln nu``, with ``lam = 1 - ef``, ``z = lam h``:
        ``nu_i^lam h (phi1(z) - phi2(z))`` to the left point and ``nu_i^lam h
        phi2(z)`` to the right, ``phi1 = (e^z - 1)/z``, ``phi2 = int_0^1 t
        e^{z t} dt``, by series near ``z = 0``.  `ef = 1` is the trapezoid in
        `ln nu`; `ef = 0` is a density sampled at the points.

        ⚠ AND RICHARDSON ON TOP, on an even number of intervals: ``(4 Q_h -
        Q_2h) / 3`` with ``Q_2h`` on every other point.  The product rule is
        exact for a flat response but carries ``h^2 g''`` where the response
        BENDS, and at `ef != 1` that term does not telescope to the band
        ends (it does at `ef = 1`): measured on an RC under 1/f^0.8, +6.85e-5
        / +1.71e-5 / +4.28e-6 / +1.07e-6 at 20 / 40 / 80 / 160 per decade --
        4.0x per halving, a clean `h^2`, which the combination removes.  The
        weights stay positive (Simpson's pattern at `ef = 1`), and it stays
        exact wherever the product rule is.

        History: `doc/shooting_history.md`, `PAC._coloured_covariance`."""
        nus = np.asarray(nus, dtype=float)
        if richardson and nus.size >= 3 and (nus.size - 1) % 2 == 0:
            q1 = cls._power_law_weights(nus, ef, richardson=False)
            q2 = np.zeros(nus.size)
            q2[::2] = cls._power_law_weights(nus[::2], ef, richardson=False)
            return (4.0 * q1 - q2) / 3.0
        wl, wr = cls._interval_weights(nus[:-1], nus[1:], ef)
        q = np.zeros(nus.size)
        q[:-1] += wl
        q[1:] += wr
        return q

    @staticmethod
    def _interval_weights(lo, hi, ef):
        """The product rule on intervals ``[lo, hi]`` (arrays): the weights
        of the left and right values -- see `_power_law_weights`."""
        lo = np.asarray(lo, dtype=float)
        h = np.log(np.asarray(hi, dtype=float) / lo)
        lam = 1.0 - float(ef)
        z = lam * h
        small = np.abs(z) < 1e-3
        zs = np.where(small, 1.0, z)
        phi1 = np.where(small, 1.0 + z / 2.0 + z ** 2 / 6.0 + z ** 3 / 24.0
                        + z ** 4 / 120.0, np.expm1(zs) / zs)
        phi2 = np.where(small, 0.5 + z / 3.0 + z ** 2 / 8.0 + z ** 3 / 30.0
                        + z ** 4 / 144.0, (np.exp(zs) * (zs - 1.0) + 1.0) / zs ** 2)
        base = lo ** lam * h
        return base * (phi1 - phi2), base * phi2

    @staticmethod
    def _loglog_integral(S, fs):
        """``int S df`` over the grid `fs` (last axis of `S`), `S` a power
        law between neighbouring points: on `[f1, f2]`, ``S = S1
        (f/f1)^p`` through both ends, integrated exactly --
        ``S1 f1 ln(r) (e^z - 1)/z``, ``z = ln(S2 f2 / (S1 f1))``, stable at
        ``z -> 0`` (1/f).  An interval with an end that is not positive
        (no power law passes through it) takes the linear trapezoid.

        History: `doc/shooting_history.md`, `PAC.sampled_variance`."""
        S = np.asarray(S, dtype=float)
        f1, f2 = fs[:-1], fs[1:]
        S1, S2 = S[..., :-1], S[..., 1:]
        L = np.log(f2 / f1)
        ok = (S1 > 0.0) & (S2 > 0.0)
        with np.errstate(divide='ignore', invalid='ignore'):
            z = np.log(np.where(ok, S2, 1.0) / np.where(ok, S1, 1.0)) + L
            phi = np.where(np.abs(z) < 1e-8, 1.0 + 0.5 * z,
                           np.expm1(z) / np.where(z == 0.0, 1.0, z))
        pw = S1 * f1 * L * phi
        lin = 0.5 * (S1 + S2) * (f2 - f1)
        return np.sum(np.where(ok, pw, lin), axis=-1)

    def _cached_root(self, cy_at):
        """`w -> _psd_sqrt(cy_at(w))`, cached per frequency: the sample
        series reads the same `|f + n f0|` for every instant."""
        cache = {}

        def root(w):
            k = float(w)
            if k not in cache:
                cache[k] = self._psd_sqrt(cy_at(k))
            return cache[k]
        return root

    def _perband_root_sampler(self, pss, key, states, wlo, f0, L):
        """`w -> (K, m, r)`: the COLUMNS of one per-band element at the
        injection points, for the sample series -- the element alone
        (`_one_element_cy`), cached per frequency, and classified once
        (`_perband_classify`):

          STATIONARY  the same at every point: one point, broadcast;
          SEPARABLE   ``C(x, w) = C(x, w_ref) s(w)``: the columns per point
                      at `w_ref` once, times ``sqrt(s(w))`` from one point;
          otherwise   the element at every point for every frequency.

        The columns are the element's SIGNED amplitudes where it states
        them (`_perband_amplitudes`), else the root of its PSD; `.signed`
        on the returned function says which.

        History: `doc/shooting_history.md`, `PAC._perband_root_sampler`."""
        wref = 2.0 * np.pi * f0
        wt = sorted({float(wlo), wref, 20.0 * np.pi * f0,
                     2.0 * np.pi * (float(L) + 0.5) * f0})
        kind, Cs, Ws, mode = self._perband_classify(pss, key, states, wt)
        K = len(states)
        if kind == 'stationary':
            one = [states[0]]
            root = self._cached_root(lambda w: np.broadcast_to(
                self._one_element_cy(pss, key, w, one), (K,) + Cs[0].shape[1:]))
        elif kind == 'separable':
            Cref = Cs[wt.index(wref)]
            Wref = self._psd_sqrt(Cref) if Ws is None else Ws[wt.index(wref)]
            jr, pi_, qi = np.unravel_index(int(np.argmax(np.abs(Cref))), Cref.shape)
            xr, cref = [states[jr]], complex(Cref[jr, pi_, qi])
            cache = {}

            def root(w):
                k = float(w)
                if k not in cache:
                    c = self._one_element_cy(pss, key, k, xr)[0, pi_, qi]
                    cache[k] = np.sqrt(max(float(np.real(c / cref)), 0.0)) * Wref
                return cache[k]
        else:
            root = self._perband_root(pss, key, states, mode)
        root.signed = Ws is not None
        return root

    @staticmethod
    def _psd_sqrt(Cs):
        """Symmetric PSD square roots of a stack `(..., n, n)`."""
        Cs = np.asarray(Cs, dtype=complex)
        Cs = 0.5 * (Cs + np.conj(np.swapaxes(Cs, -1, -2)))
        lam, U = np.linalg.eigh(Cs)
        return np.einsum('...ik,...k,...jk->...ij', U,
                         np.sqrt(np.clip(np.real(lam), 0.0, None)), U.conj())

    def _colour_groups(self, pss, model, states, wlo, f0, L, what,
                       warn_touch=True):
        """The coloured components of `model` as unit processes through
        their own columns: `('fixed', G, s)` -- columns `G (K, m, r)` at
        `states` and a power weight `s(nu)` (a uniform power law, its
        element's SIGNED amplitudes where stated, else the root of its PSD)
        -- or `('band', root, None)`, the columns per band frequency (a
        power law whose exponent varies across its entries; a per-band
        colour, `_perband_root_sampler`).  A component factored by the root
        of its PSD whose PSD TOUCHES ZERO along the orbit is warned on: if
        its modulation changes sign there, that root is the `|m|` process
        (`warn_touch`; the sample series has never warned it).  Every band
        root is cached per frequency (`_cached_root`): the sample series
        reads the same `|f + n f0|` at every instant.  The FIXED groups
        before the BAND ones is the order the sample series sums them in."""
        signed = getattr(model, 'amplitude', None) or {}
        self._warn_signed_unused(model, 'PAC.%s' % what)

        def touches(C):
            ## the necessary condition for a sign change, as the pnoise fold
            ## asks it: a diagonal entry that falls to 1e-2 of its maximum
            d = np.abs(np.real(np.diagonal(np.asarray(C), axis1=-2, axis2=-1)))
            dmax = d.max(axis=0)
            return bool(np.any((dmax > 0) & (d.min(axis=0) <= 1e-2 * dmax)))
        groups, blind = [], []
        for key, Bc, EF in model.flicker:
            ef = self._uniform_exponent(Bc, EF)
            W = signed.get(key)
            if ef is not None:
                if W is None:
                    if warn_touch and touches(Bc):
                        blind.append(key)
                    W = self._psd_sqrt(Bc)
                groups.append(('fixed', np.asarray(W, dtype=complex),
                               lambda nu, ef=ef, w1=model.w1:
                               (w1 / np.asarray(nu, dtype=float)) ** ef))
            else:
                ## (signed columns grouped by their own exponents, where
                ## they carry one each: a uniform power law per group)
                split = (self._exponent_columns(Bc, EF, W) if W is not None
                         else None)
                if split is not None:
                    for Wg, efg in split:
                        groups.append(('fixed', Wg,
                                       lambda nu, ef=efg, w1=model.w1:
                                       (w1 / np.asarray(nu, dtype=float)) ** ef))
                    continue
                if warn_touch and touches(Bc):
                    blind.append(key)
                groups.append(('band', self._cached_root(
                    lambda w, Bc=Bc, EF=EF, w1=model.w1: Bc * (w1 / w) ** EF),
                    None))
        for key in model.perband:
            root = self._perband_root_sampler(pss, key, states, wlo, f0, L)
            ## (the test costs a CY evaluation per point: only when it can
            ## warn, and an element's signed amplitudes carry the sign)
            if warn_touch and not root.signed and touches(self._one_element_cy(
                    pss, key, 2.0 * np.pi * f0, states)):
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

    def _refuse_driven(self, pss, what):
        if not getattr(pss, 'autonomous', False):
            raise ValueError(
                'PAC.%s: phase diffusion is a property of a '
                "FREE-RUNNING oscillator. A driven circuit's phase is its "
                "source's, and its noise is pnoise's problem, not this one."
                % what)

    @staticmethod
    def _period_weights(tms, nsamp, T, pss=None):
        """Periodic TRAPEZOID weights for samples at `tms[0..nsamp-1]` over
        a period `T`: `w_j = (g_j + g_{j-1}) / 2` with `g_j` the gap to the
        next sample and the last gap closing the period.

        ⚠ NOT `h = diff(times)`: that is the LEFT RECTANGLE rule, spectrally
        accurate on a uniform periodic grid but FIRST order on a NON-UNIFORM
        one (its error is (1/2) integral h'(t) y(t) dt, not zero).  On a
        uniform grid `0.5 h + 0.5 h == h` exactly.  ⚠ AND THE TRAPEZOID IS
        ITSELF A SECOND-ORDER CAP on a smoothly varying grid: with `pss`
        given and a non-uniform grid, these are the periodic cubic-spline
        weights of `periodic_spline_weights`, broken at the nodes of landed
        events; a uniform grid is unchanged.

        History: `doc/shooting_history.md`, `PAC._period_weights`."""
        tms = np.asarray(tms, dtype=float).ravel()
        n = int(nsamp)
        t = tms[:n]
        g = np.empty(n)
        g[:-1] = t[1:] - t[:-1]
        g[-1] = float(T) + t[0] - t[-1]
        if (pss is not None and n >= 4
                and float(np.max(g)) / float(np.min(g)) - 1.0 > pss.UNIFORM_GRID_TOL):
            return periodic_spline_weights(t, T, pss._event_nodes(t, T))
        w = 0.5 * g
        w[1:] += 0.5 * g[:-1]
        w[0] += 0.5 * g[-1]
        return w

    def _ppv_states(self, pss):
        """The full-width orbit states the PPV samples and the Floquet modes
        live on, one per sample: the monodromy twin's orbit where a twin
        serves them (trap, euler), else the solve's own."""
        tw = pss.monodromy_twin() if hasattr(pss, 'monodromy_twin') else pss
        return self._orbit_states(tw)

    def _modulated_present(self, pss):
        """Whether a noise source of the circuit follows the orbit -- the
        bias dependence `_cy_reduced` refuses, asked without refusing."""
        try:
            self._cy_reduced(pss, 2.0 * np.pi / float(pss.period))
        except NotImplementedError:
            return True
        return False

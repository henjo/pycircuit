"""A converged period kept factored (`FactoredPeriod`, one class per kind), and
what one walk of the period produced (`_PeriodWalk`).
"""
import numpy as np
from ._steps import _GLMStep
from ._steps import _LMMStep


class FactoredPeriod(object):
    """One converged period, kept FACTORED -- the hook PAC/PPV/pnoise share.

    `PSS.solve` throws its factorisations away: they live inside the Newton's
    `build` closure and go out of scope with it.  Every small-signal analysis
    over a periodic operating point needs the same thing back -- products
    with the monodromy, never the monodromy itself -- so this is what
    `PSS.factored_period()` hands out.

    ⚠ `matvec` IS THE WHOLE INTERFACE, and deliberately so.  The withdrawn
    `PAC` read `pss.Jtvec` / `pss.Cvec` and rebuilt the `(NM)x(NM)` operator
    from them, which is where its 419.5 GiB went.  Two things are wrong with
    that route and only one of them is the memory: those two lists are
    written by `_traverse` and `_traverse_solved_history` and by NEITHER
    factored traversal, so after `solve(matrix_free=True)` they are stale or
    absent -- an analysis reading them would rebuild an operator for a
    DIFFERENT trajectory than the one that converged, silently.

    ⚠ AND THE REBUILT OPERATOR WAS EULER-SHAPED.  See
    `test_the_pac_L_is_backward_euler_only`: the old `L` has two terms per
    row, so for `trap` or `gear` it is not the discretisation the trajectory
    was produced by and `L^-1 B` is not the monodromy at all -- measured,
    spectral radius 0 against the analytic 0.8546.  A `matvec` cannot make
    that mistake, because each step carries its own `(alphas, b)`.
    """

    __slots__ = ('kind', 'opening', 'steps', 'x_last', 'x_prev', 'width',
                 'times', 'T', 'open_at_x0', '_pss')
    ## 'glm' is the MULTIVALUE kind: width r*m, see `factored_period_glm`.

    ## ONE CLASS PER KIND (2026-09-24): `FactoredPeriod(kind, ...)` builds the
    ## subclass that kind names (`_PERIOD_KINDS`), and what differs between
    ## the kinds -- how a direction seeds the per-step state, what the map
    ## reads out of it, where a costate injection lands -- is that class's
    ## own methods, where it was a four-way branch in each of seven.
    ## Consumers ask what a map IS through these flags, not its kind string.
    is_plain = is_pair = is_stage = is_glm = False

    def __new__(cls, kind=None, *args, **kwargs):
        if cls is FactoredPeriod and kind is not None:
            cls = _PERIOD_KINDS[kind]
        return object.__new__(cls)

    def __init__(self, kind, opening, steps, x_last, x_prev, pss,
                 times=None, T=None, open_at_x0=False):
        self.kind, self.opening, self.steps = kind, opening, steps
        self.x_last, self.x_prev, self._pss = x_last, x_prev, pss
        ## the grid the steps were taken on -- a forced replay needs the
        ## time of each step to evaluate `exp(j w t)` there, and reading it
        ## off `pss.times` later is exactly the parallel-indexing trap the
        ## final replay's `(t, h)` pairing was rewritten to remove
        self.times, self.T = times, T
        ## whether the period opened AT `x_0` -- no manufacturing step.  A
        ## driven replay needs to know: the manufacturing step is not in
        ## `steps`, so its share of the source is not applied, and that
        ## costs an order.  See `PAC.solve`.
        self.open_at_x0 = bool(open_at_x0)
        m = pss.cir.n - 1
        ## the solved-history map acts on the PAIR, the plain map on one
        ## state -- the same distinction `_monodromy` carries
        self.width = 2 * m if self.is_pair else m

    def matvec(self, v):
        """`M v`, real or complex, replaying the stored factors (the one
        replay every kind shares: `PSS._replay`)."""
        return self._pss._replay(self, v)

    def matvec_transposed(self, v, collect=False, inject=None):
        """`M^T v` -- see `_monodromy_matvec_transposed{,_plain}`.

        ⚠ B8: the PLAIN path now has one too, for the ONE-STEP companions.
        It used to refuse outright, which made every adjoint surface --
        `ppv`, PAC, `pnoise`, `covariance` -- Gear-2 only.

        ⚠ `collect` AND `inject` FORWARD TO BOTH, which is what lets the
        callers stop naming a map. They used to reach past this method to
        `_monodromy_matvec_transposed` DIRECTLY, so B8 shipped without
        reaching any of them -- the machinery existed and every surface
        still refused. A caller that goes through here gets whichever
        recursion its `kind` calls for and cannot acquire a Gear-2
        assumption by accident.

        ⚠ THE COLLECTED STATE HAS THE MAP'S OWN WIDTH: `2m` for the pair,
        `m` for the plain one. `st[:m]` is the differential block under
        both, which is the slice every consumer wants.
        """
        return self._pss._replay_transposed(self, v, collect=collect,
                                            inject=inject)

    ## -- what differs between the kinds: one subclass each -----------------
    ##
    ## (2026-09-23) Every kind's map is `extract . step_N ... step_1 . seed`:
    ## the steps are objects with one algebra (`_StageStep`, `_LMMStep`,
    ## `_GLMStep`: `solve`, `adjoint`, `sources`, `source_adjoint`), and what
    ## is left per kind is how a direction seeds the per-step state, what the
    ## map reads out of it, and where a costate injection lands.  The replays
    ## (`PSS._replay`, `_replay_transposed`, `_forced_replay`,
    ## `_forced_replay_transposed`, `_sideband_forced`) are one function each.
    ## The interface, defined by `_PlainPeriod`, `_PairPeriod`,
    ## `_StagePeriod` and `_GLMPeriod`:

    def step_objects(self):
        """The steps as objects that know their own algebra."""
        raise NotImplementedError

    def seed(self, v):
        """The per-step state the map starts from, for a direction `v` of
        width `self.width`."""
        raise NotImplementedError

    def extract(self, c):
        """What the map reads out of the final per-step state (width
        `self.width`)."""
        raise NotImplementedError

    def node(self, c):
        """The circuit state at the node a per-step state ends on (width
        `m`) -- what a forced replay collects."""
        raise NotImplementedError

    def extract_T(self, v):
        """The adjoint of `extract`: the final adjoint state for a seed
        `v`."""
        raise NotImplementedError

    def seed_T(self, w):
        """The adjoint of `seed`: back to a direction of width
        `self.width`."""
        raise NotImplementedError

    def inject(self, w, x):
        """A costate injection `x` at a node."""
        raise NotImplementedError

    def collected(self, w):
        """The adjoint state at a node, as `matvec_transposed(collect=True)`
        returns it; ``st[:m]`` is the circuit block under every map."""
        raise NotImplementedError


class _LMMPeriod(FactoredPeriod):
    """A linear-multistep map, plain or gear's pair.  Its steps are stored
    records ``(lu, C_new, alphas, b)`` and become `_LMMStep`s with the
    capacitance ring the forward pass saw -- rebuilt here, `N` references to
    matrices that already exist -- and each step's end time.  The per-step
    state is ``(P_n, P_{n-1}, Pq_n)`` for both; they differ in how a
    direction seeds it and what the map reads out.  A costate injection
    lands on the circuit state ``P_n``."""
    __slots__ = ()

    def step_objects(self):
        ring = self._ring()
        tms = self.times
        out = []
        for j, st in enumerate(self.steps):
            out.append(_LMMStep(st, ring[0], ring[1],
                                None if tms is None else tms[j + 1]))
            ring = [st[1], ring[0]]
        return out

    def node(self, c):
        return c[0]

    def inject(self, w, x):
        return (w[0] + x, w[1], w[2])


class _PlainPeriod(_LMMPeriod):
    """The PLAIN map ('plain', width `m`): one entering state, both ring
    slots opened on it.

    Its `seed` takes ``Pq`` from THE OPENING PAIR, not the loop's -- the
    walk opens it right after the MANUFACTURING step, order-dropped to Euler
    (``b = 0``), so in practice at zero; using the loop's made the map 100 %
    wrong for `trap`.  Plus the consistent-``iq_0`` seed (`_pq_seed_at_x0`,
    `theta`), `None` for every other method."""
    __slots__ = ()
    is_plain = True

    def _ring(self):
        return [self.opening[0], self.opening[0]]

    def seed(self, v):
        m = self._pss.cir.n - 1
        C_open, a_open, b_open, pq_open = self.opening
        Pq = (a_open[0] * (C_open @ v) if b_open
              else np.zeros(m, dtype=v.dtype))
        if pq_open is not None:
            Pq = Pq + pq_open @ v
        return (v.copy(), v.copy(), Pq)

    def extract(self, c):
        return c[0]

    def extract_T(self, v):
        m = self._pss.cir.n - 1
        return (v.copy(), np.zeros(m, dtype=v.dtype),
                np.zeros(m, dtype=v.dtype))

    def seed_T(self, w):
        ## ⚠ THE SEED'S COMPANION CURRENT IS ONLY DISCARDABLE WHEN IT DOES
        ## NOT DEPEND ON THE SEED: the forward map opens at `P = v`, `Pq =
        ## (a_0 C_open [b_open] + pq_open) v`, so the transpose closes at
        ## `w1 + w2 + (...)^T w3` (`w2`, the adjoint of the carried slot, is
        ## zero on a one-step map)
        C_open, a_open, b_open, pq_open = self.opening
        out = w[0] + w[1]
        if b_open:
            out = out + a_open[0] * (np.asarray(C_open).T @ w[2])
        if pq_open is not None:
            out = out + pq_open.T @ w[2]
        return out

    def collected(self, w):
        ## the circuit block (`m`): the plain map's `Pq` adjoint is a
        ## companion term, not a state
        return w[0].copy()


class _PairPeriod(_LMMPeriod):
    """Gear's solved-history PAIR ('solved_history', width `2m`): the
    direction is ``(v_0, v_{-1})``, the companion current opens at zero
    (``b = 0``), and the map reads out -- and collects -- both blocks.
    ``st[:m]`` is the circuit block, as under every map."""
    __slots__ = ()
    is_pair = True

    def _ring(self):
        return list(self.opening)

    def seed(self, v):
        m = self._pss.cir.n - 1
        return (v[:m].copy(), v[m:].copy(), np.zeros(m, dtype=v.dtype))

    def extract(self, c):
        return np.concatenate((c[0], c[1]))

    extract_T = seed

    def seed_T(self, w):
        return np.concatenate((w[0], w[1]))

    def collected(self, w):
        return np.concatenate((w[0].copy(), w[1].copy()))


class _StagePeriod(FactoredPeriod):
    """A Runge-Kutta stage map ('full' coupled, 'dirk' lower triangular,
    width `m`): self-starting, stored as `_StageStep`s, and the per-step
    state IS the direction -- nothing to seed, read out or split."""
    __slots__ = ()
    is_stage = True

    def step_objects(self):
        return self.steps

    def seed(self, v):
        return v

    def extract(self, c):
        return c

    def node(self, c):
        return c

    def extract_T(self, v):
        return v

    def seed_T(self, w):
        return w

    def inject(self, w, x):
        return w + x

    def collected(self, w):
        return w.copy()


class _GLMPeriod(FactoredPeriod):
    """A Nordsieck GLM's map ('glm', width ``r*m``): the per-step state is
    the Nordsieck blocks as columns; a costate injection lands on the whole
    Nordsieck vector; the circuit state at a node is its first block."""
    __slots__ = ()
    is_glm = True

    def step_objects(self):
        return [_GLMStep(st) for st in self.steps]

    def seed(self, v):
        m = self._pss.cir.n - 1
        return [v[k * m:(k + 1) * m].reshape(m, 1)
                for k in range(v.shape[0] // m)]

    def extract(self, c):
        return np.concatenate([np.asarray(Pk).ravel() for Pk in c])

    def node(self, c):
        return c[0]

    def extract_T(self, v):
        m = self._pss.cir.n - 1
        return [v[k * m:(k + 1) * m].copy() for k in range(v.shape[0] // m)]

    def seed_T(self, w):
        return np.concatenate(w)

    def inject(self, w, x):
        m = self._pss.cir.n - 1
        x = np.asarray(x).ravel()
        return [w[k] + x[k * m:(k + 1) * m] for k in range(len(w))]

    def collected(self, w):
        return np.concatenate([x.copy() for x in w])


_PERIOD_KINDS = {'plain': _PlainPeriod, 'solved_history': _PairPeriod,
                 'full': _StagePeriod, 'dirk': _StagePeriod,
                 'glm': _GLMPeriod}


class _PeriodWalk(object):
    """What one walk of the period produced (`PSS._walk_lmm`,
    `PSS._walk_stage`): the opened and final states (`x_prev` the one a step
    before the end), the dense monodromy `P` (a two-entry ring on the
    multistep walk), the period column `Pt`, the event columns `Pk`, the kept
    steps and what a replay seeds from (`opening`).  A field the walk was not
    asked for is None.

    `kind` is the map's (`PSS._map_kind`), and the methods read the walk in
    the MAP'S coordinates -- gear's pair stacks ``(x_0, x_{-1})`` and ends on
    ``(x_{N-1}, x_{N-2})`` -- so a caller needs no branch per kind.
    `fp_kind` is the `FactoredPeriod` kind the kept steps make."""
    __slots__ = ('kind', 'fp_kind', 'z', 'x0', 'x_end', 'x_prev', 'P', 'Pt',
                 'Pk', 'steps', 'opening', 'open_at_x0', 'width')

    def __init__(self, kind=None, fp_kind=None, z=None, x0=None, x_end=None,
                 x_prev=None, P=None, Pt=None, Pk=None, steps=None,
                 opening=None, open_at_x0=False, width=None):
        self.kind, self.fp_kind, self.z = kind, fp_kind, z
        self.x0, self.x_end, self.x_prev = x0, x_end, x_prev
        self.P, self.Pt, self.Pk = P, Pt, Pk
        self.steps, self.opening = steps, opening
        self.open_at_x0, self.width = open_at_x0, width

    def z0(self):
        """Where the map starts: the pair's unknown itself, else the state
        the period opened at (after the manufacturing step, if any)."""
        return self.z if self.kind == 'pair' else self.x0

    def end(self):
        """Where the map ends, in the unknown's coordinates."""
        if self.kind == 'pair':
            return np.concatenate((np.asarray(self.x_end),
                                   np.asarray(self.x_prev)))
        return self.x_end

    def monodromy(self):
        """The dense map `M = d end / d z0`."""
        if self.kind == 'pair':
            return np.vstack((self.P[0], self.P[1]))
        return self.P[0] if self.kind == 'plain' else self.P

    def period_column(self):
        """`d end / dT`, one column in the map's coordinates."""
        if self.kind == 'pair':
            return np.concatenate((np.asarray(self.Pt[0]).ravel(),
                                   np.asarray(self.Pt[1]).ravel()))
        return self.Pt[0] if self.kind == 'plain' else self.Pt

    def factored(self, pss, times=None, T=None):
        """The kept steps as a `FactoredPeriod` (a GLM's acts on the
        Nordsieck state, `width` wide)."""
        fp = FactoredPeriod(self.fp_kind, self.opening, self.steps,
                            self.x_end,
                            (self.x0 if self.kind in ('plain', 'glm')
                             else self.x_prev),
                            pss, times=times, T=T, open_at_x0=self.open_at_x0)
        if self.width is not None:
            fp.width = self.width
        return fp

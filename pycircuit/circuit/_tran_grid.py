"""The transient's uniform output grid: `resample_uniform`, the interpolant,
and `OutputGrid`, which makes `outputstep`'s grid as the run goes, so a run
keeps the grid and its last few points rather than every accepted state.  A
theme of `Transient` (see `transient.py`, which re-exports
`resample_uniform`).
"""

import numpy as np

#: `OutputGrid` makes a block of the grid once the run's finished points hold
#: this many bytes; a run below it is one block, made at its end.
GRID_BLOCK_BYTES = 1 << 26
#: A point is finished, never written again, once this many follow it: the
#: periodic shifts (`_apply_periodic_shifts`) rewrite the newest three.
FINISHED_AFTER = 3


def _uniform_grid(t0, t1, step):
    """`step`'s grid from `t0` to `t1`, `t1` included when the step divides
    the interval, nothing beyond it.

    `t1 + step/2` so a grid that divides the interval exactly still includes
    the endpoint, without inventing a point beyond it."""
    grid = np.arange(t0, t1 + float(step) / 2.0, float(step))
    return grid[grid <= t1 * (1 + 1e-12)]


def resample_uniform(t, x, npoints=None, step=None, grid=None):
    """Interpolate a transient result onto a UNIFORM time grid.

    A transient returns the solver's own adaptive points, so their spacing is
    whatever the step controller chose.  Anything that needs a uniform grid,
    above all an FFT, resamples with this rather than with `np.interp`, which
    is LINEAR and throws away an order of accuracy before the transform.

    QUADRATIC, to match the integrator.  Three-point Lagrange through the
    interval's neighbours, the same reasoning `TLine._interpolate_history` gives
    for its own history lookup: a first-order interpolant feeding a second-order
    method injects error the method never made (measured 14-41x below linear
    where the grid resolves the signal, 1.3x where it barely does).  Where the
    adaptive grid is barely resolving the signal, no interpolant recovers what
    was not sampled, and the right fix there is a smaller `max_step`, not a
    cleverer resample.

    Give exactly one of ``npoints``, ``step`` or ``grid``.

    ``grid`` takes an explicit set of output times, for a caller whose window
    this function cannot otherwise express -- an FFT that needs a TRAILING
    window with ``endpoint=False``, because dropping the duplicated period
    boundary is what keeps the tones exactly on bins.  Neither ``npoints``
    (which spans ``t[0]..t[-1]`` inclusive) nor ``step`` (which starts at
    ``t[0]``) can say that.

    History: `doc/transient_history.md`, `resample_uniform`.
    """
    numpy = np

    t = numpy.asarray(t, dtype=float)
    x = numpy.asarray(x)
    given = sum(spec is not None for spec in (npoints, step, grid))
    if given != 1:
        raise ValueError('give exactly one of npoints, step or grid')
    if t.size < 3:
        raise ValueError('need at least 3 points to interpolate quadratically, '
                         f'got {t.size}')

    if grid is not None:
        grid = numpy.asarray(grid, dtype=float)
    elif npoints is not None:
        grid = numpy.linspace(t[0], t[-1], int(npoints))
    else:
        grid = _uniform_grid(t[0], t[-1], step)

    ## The interval containing each output point, clamped so the three-point
    ## stencil stays inside the data at both ends.
    idx = numpy.clip(numpy.searchsorted(t, grid) - 1, 1, t.size - 2)
    t0, t1, t2 = t[idx - 1], t[idx], t[idx + 1]
    L0 = (grid - t1) * (grid - t2) / ((t0 - t1) * (t0 - t2))
    L1 = (grid - t0) * (grid - t2) / ((t1 - t0) * (t1 - t2))
    L2 = (grid - t0) * (grid - t1) / ((t2 - t0) * (t2 - t1))

    if x.ndim == 1:
        out = x[idx - 1] * L0 + x[idx] * L1 + x[idx + 1] * L2
    else:
        ## Rows are unknowns, columns are time -- the shape a CircuitResult holds.
        out = (x[:, idx - 1] * L0 + x[:, idx] * L1 + x[:, idx + 1] * L2)
    return grid, out


class OutputGrid:
    """`outputstep`'s grid made as the run goes: `resample_uniform(t, X,
    step=outputstep)` of the whole run, bit for bit, without the whole run
    held.

    The run appends each accepted state to `X` and its time to `T` (the
    times after `X[0]`'s, which is at 0) and calls `advance`.  A point is
    FINISHED once `FINISHED_AFTER` more follow it: the periodic shifts
    (`_apply_periodic_shifts`) rewrite the newest three in place, and nothing
    else writes a past point.  When the finished points not yet made into the
    grid hold `GRID_BLOCK_BYTES`, the grid points at or before the last of
    them are made by `resample_uniform` on those points alone, and all but the
    last two are dropped from `X` (their entries None: `len(X)` and the live
    entries' indices hold; every reader of `X` takes its last three).

    Bit for bit, because each grid point gets the run's own stencil and the
    same expressions on it.  The stencil of a grid point in `(t[j-1], t[j]]`
    is points `j-2, j-1, j`; a block from point `f` to `l` makes the grid
    points in `(t[f+1], t[l]]` (from `t[0]` in the first), whose stencils lie
    inside it, its search finding the same interval and its clamp the run's
    at the run's two ends.  The grid's times are numpy's `arange` from 0 by
    the step, which a longer `arange` starts with (checked at `finish`).
    """

    def __init__(self, X, T, step, tend):
        step = float(step)
        tend = np.float64(tend)
        self.X, self.T, self.step = X, T, step
        ## the grid past `tend`, so a run ending a rounding beyond it is
        ## inside; the grid of a run that ends at `tend`, its size
        self.grid = np.arange(np.float64(0.0), tend + step / 2.0 + 3.0 * step, step)
        self.size = len(_uniform_grid(np.float64(0.0), tend, step))
        self.n, self.row = X[0].size, X[0].nbytes
        self.first = 0
        self.made = 0
        self.out = None

    def advance(self):
        """After a point is appended to `X` and `T`: a block made when the
        finished points fill one."""
        last = len(self.X) - 1 - FINISHED_AFTER
        if last - self.first >= 2 and (last - self.first + 1) * self.row >= GRID_BLOCK_BYTES:
            self._block(last)

    def _times(self, last):
        """The times of points `first..last`."""
        T, first = self.T, self.first
        return np.array([0.0] + T[:last] if first == 0 else T[first - 1:last])

    def _block(self, last):
        """The grid points at or before `t[last]` made from points
        `first..last`, all but the last two of them dropped."""
        X, first = self.X, self.first
        t = self._times(last)
        stop = int(np.searchsorted(self.grid, t[-1], side='right'))
        if stop > self.made:
            if self.out is None:
                self.out = np.empty((self.n, self.size))
            _g, y = resample_uniform(t, np.array(X[first:last + 1]).T,
                                     grid=self.grid[self.made:stop])
            self.out[:, self.made:stop] = y
            self.made = stop
        for k in range(first, last - 1):
            X[k] = None
        self.first = last - 1

    def finish(self):
        """`(grid, values)`: the rest of the grid made from the points left,
        as `resample_uniform` makes the whole run's."""
        X, T = self.X, self.T
        if len(X) < 3:
            ## (its refusal, as the whole run's)
            return resample_uniform(np.array([0.0] + T), np.array(X).T, step=self.step)
        grid = _uniform_grid(np.float64(0.0), np.float64(T[-1]), self.step)
        if len(grid) > len(self.grid) or not np.array_equal(grid, self.grid[:len(grid)]):
            raise RuntimeError('outputstep: the grid made during the run is not the run\'s')
        if self.made > len(grid):
            raise RuntimeError('outputstep: the grid made during the run passed its end')
        last = len(X) - 1
        t = self._times(last)
        _g, y = resample_uniform(t, np.array(X[self.first:]).T, grid=grid[self.made:])
        if self.out is None:
            return grid, y
        out = self.out
        if out.shape[1] != len(grid):
            out = np.empty((self.n, len(grid)))
            out[:, :self.made] = self.out[:, :self.made]
        out[:, self.made:] = y
        for k in range(self.first, last + 1):
            X[k] = None
        self.out = None
        return grid, out

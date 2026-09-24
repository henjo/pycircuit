"""The shooting analyses: the periodic steady state (`PSS`), the periodic
small-signal and noise analyses over it (`PAC`), and probe shooting
(`ProbeShooting`).

Until 2026-09-24 this was one module, `shooting.py` (18,985 lines).  It is a
package of:

    pss.py            PSS: the class, its parameters, `solve` in its phases
    _pss_newton.py    its outer solves: free-period and matrix-free Newton
    _pss_grids.py     its grids: event grid, LTE grid and fold, period grid
    _pss_events.py    its state events as Newton unknowns
    _pss_inner.py     its inner transient and the circuit at a point
    _pss_walks.py     its period walks, dense and factored
    _pss_replays.py   its factored period and replays
    _pss_ppv.py       its adjoint side: PPV, Floquet modes
    _pss_accuracy.py  its accuracy checks: twin, grid_error, warping estimate
    _pss_periodic.py  its states defined up to a modulus (idtmod)
    pac.py            PAC and SidebandResponse
    probe.py          ProbeShooting
    events.py         EventColumns
    diagnostics.py    topological_index, algebraic_conditioning, ...
    _factored.py      FactoredPeriod (one class per kind), _PeriodWalk
    _steps.py         the step objects and their recursions
    _numerics.py      small numerical helpers

This file exports the package's public names; a private helper is imported
from the module that defines it.
"""
from .pss import AUTONOMOUS_U_TOL, PSS
from .pac import PAC, SidebandResponse
from .probe import ProbeShooting
from .events import EventColumns
from ._factored import FactoredPeriod
from .diagnostics import (algebraic_conditioning, noise_enters_constraints,
                          topological_index)
from ._numerics import freq_analysis, periodic_spline_weights

__all__ = ['AUTONOMOUS_U_TOL', 'EventColumns', 'FactoredPeriod', 'PAC', 'PSS',
           'ProbeShooting', 'SidebandResponse', 'algebraic_conditioning',
           'freq_analysis', 'noise_enters_constraints',
           'periodic_spline_weights', 'topological_index']

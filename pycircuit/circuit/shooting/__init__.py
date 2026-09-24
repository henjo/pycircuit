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

Every name the single module had is here too, so `shooting.X` and
`from pycircuit.circuit.shooting import *` work as they did (`__all__` is
that module's public names, exactly).
"""
## ruff: noqa: F401, F403, F405 -- this file IS the re-export: every name the
## single module had, kept importable from here
## the single module's own imports, re-exported as its namespace had them
from pycircuit.post import InternalResultDict
from pycircuit.circuit.circuit import gnd
from pycircuit.circuit.analysis import *
from copy import copy
import warnings
import pycircuit.circuit.analysis as analysis
import numpy as np

from ._numerics import (freq_analysis, _complex_solve, _complex_solve_transposed,
                         _arnoldi_gmres, periodic_spline_weights,
                         _cx_collect, _lu_solve_split, _sla_lu_solve)
from ._steps import (_StageStep, _butcher, _lmm_recursion, _LMMStep, _GLMStep, _glm_step)
from ._factored import (FactoredPeriod, _LMMPeriod, _PlainPeriod, _PairPeriod, _StagePeriod,
                         _GLMPeriod, _PERIOD_KINDS, _PeriodWalk)
from .diagnostics import (noise_enters_constraints, _TI_CAPACITIVE, _TI_VOLTAGE, _TI_INDUCTIVE,
                           _TI_CURRENT, _TI_RESISTIVE,
                           topological_index, _limit_state_snapshot,
                           _limit_state_restore,
                           algebraic_conditioning)
from .events import (EventColumns)
from .pss import (AUTONOMOUS_U_TOL, _SolveRun, PSS)
from .pac import (SidebandResponse, _output_weights, _output_row, PAC)
from .probe import (ProbeShooting)
from ._pss_newton import _ShootingNewton
from ._pss_grids import _PeriodGrids
from ._pss_events import _StateEvents
from ._pss_inner import _InnerTransient
from ._pss_walks import _PeriodWalks
from ._pss_replays import _FactoredReplays
from ._pss_ppv import _PPVFloquet
from ._pss_accuracy import _AccuracyChecks
from ._pss_periodic import _PeriodicStates

__all__ = ['AUTONOMOUS_U_TOL', 'Analysis', 'C', 'Circuit', 'CircuitResult', 'Diode',
           'EventColumns', 'FactoredPeriod', 'IS', 'IVResultDict',
           'InternalResultDict', 'L', 'NoConvergenceError', 'PAC', 'PSS',
           'Parameter', 'ParameterDict', 'ProbeShooting', 'R',
           'SidebandResponse', 'SingularMatrix', 'SubCircuit', 'VS',
           'Waveform', 'algebraic_conditioning', 'analysis',
           'analysis_kind', 'circuit', 'contextlib', 'copy', 'defaultepar',
           'freq_analysis', 'fsolve', 'gnd', 'instjoin', 'isiterable',
           'newton_tolerance_vectors', 'noise_enters_constraints', 'np',
           'numeric', 'numpy', 'periodic_spline_weights',
           'reduced_row_names', 'remove_row_col', 'sim', 'symbolic',
           'topological_index', 'types', 'warnings']

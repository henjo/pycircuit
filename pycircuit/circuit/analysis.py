# -*- coding: latin-1 -*-
# Copyright (c) 2008 Pycircuit Development Team
# See LICENSE for details.

import numpy
from pycircuit import sim
from pycircuit.utilities import Parameter, ParameterDict, isiterable
from pycircuit.circuit import Circuit, SubCircuit, VS,IS,R,C,L,Diode, gnd, \
    defaultepar, instjoin, circuit
from . import circuit
from .toolkit import symbolic
from pycircuit.post.waveform import Waveform
from pycircuit.post.result import IVResultDict
from pycircuit.post.internalresult import InternalResultDict
from copy import copy
from .toolkit import numeric
import types

class NoConvergenceError(Exception):
    pass

class SingularMatrix(Exception):
    pass


import contextlib

@contextlib.contextmanager
def analysis_kind(epar, kind):
    """Scope ``epar.analysis_kind = kind`` for the duration of a solve.

    Elements whose stamps differ per analysis (the ``Idt``/``Idtmod`` DC ic
    pin, idtmod.md sec. 5.1) read this flag from the epar their ``G``/``i``/
    ``u`` receive.  The ``finally`` restore is load-bearing: analyses share
    epar objects (``Transient`` forwards its own into its inner ``DC``; the
    JAX assemblies read ``defaultepar``), so a leaked flag pins every
    integrator for whatever runs next -- VACASK shipped exactly that failure
    through OSDI's EnableIntegration.  Nesting restores the previous value,
    not None.
    """
    if 'analysis_kind' not in epar:
        epar.append(Parameter(name='analysis_kind',
                              desc='Which analysis is evaluating the '
                                   'circuit; set transiently by analyses',
                              default=None))
    previous = epar.analysis_kind
    epar.set(analysis_kind=kind)
    try:
        yield
    finally:
        epar.set(analysis_kind=previous)

def reduced_row_names(cir, irefnode):
    """Names for the solver's rows, with the reference row removed.

    STAGE 6.  The solver works on a system with the reference node eliminated, so
    a row index it reports is a *reduced* index and cannot be looked up in the
    circuit directly.  This produces the matching name list once, so a diagnostic
    can say "'net3'" instead of "row 4".

    Returns None if the names cannot be built, which keeps every diagnostic that
    uses it optional -- a solver on a system with no circuit behind it (the
    nrsolver unit tests do exactly that) still works, it just reports indices.
    """
    try:
        names = [str(n.name) for n in cir.nodes]
        names += ['branch %s' % getattr(b, 'name', i)
                  for i, b in enumerate(cir.branches)]
    except Exception:
        return None
    if 0 <= irefnode < len(names):
        del names[irefnode]
    return names


class CircuitResult(IVResultDict, InternalResultDict):
    """Result class for analyses that returns voltages and currents"""
    def __init__(self, circuit, x, xdot = None, 
                 sweep_values=[], sweep_label='', sweep_unit=''):
        super().__init__()

        nodes = circuit.nodes

        self.circuit = circuit
        self.x = x
        self.xdot = xdot
        self.sweep_values = sweep_values
        self.sweep_label = sweep_label
        self.sweep_unit = sweep_unit

    def build_waveform(self, result, ylabel, yunit):
        if hasattr(result, '__iter__'):
            return Waveform(self.sweep_values, result,
                            ylabel = ylabel, yunit = 'V', 
                            xlabels = (self.sweep_label,), 
                            xunits=(self.sweep_unit,))
        else:
            return result

    def v(self, plus, minus=None):
        result = self.circuit.extract_v(self.x, plus, minus)

        if minus is not None:
            ylabel = 'v(%s,%s)'%(str(plus), str(minus))
        else:
            ylabel = 'v(%s)'%(str(plus))

        return self.build_waveform(result, ylabel, 'V')

    def i(self, term):
        """Return terminal current i(term)"""
        result = self.circuit.extract_i(self.x, term, xdot = self.xdot)    
        return self.build_waveform(result, 'i(%s)'%(str(term)), 'A')

def _reduce_ndarray(A, n):
    """Drop index ``n`` from every axis of a numpy array, in one pass per array.

    STAGE 7a.  The generic path below calls ``toolkit.delete`` once per axis, so a
    square Jacobian is copied TWICE -- an ``(n-1, n)`` intermediate and then the
    ``(n-1, n-1)`` result.  Copying the four surviving blocks straight into one
    output does the same work once, measured at **2.2x to 7.4x** faster from
    n=200 to n=3200 with the allocation included.

    The allocation is included deliberately: a buffer cached across calls would be
    faster still and is not safe here, because callers hold on to what they are
    given -- the transient keeps ``J`` for the step controller's LTE solve after
    the Newton step has used it, so a second call would overwrite the first
    caller's matrix.

    Why this matters at all, given it is 0.2% of a transient at n=302: it is the
    worst-SCALING component of the step, measured at n^2.53 against assembly's
    n^1.26 and dense LU's n^1.83.  At n=20000 the old form costs 12.4 s per Newton
    iteration -- more than the LU it feeds.
    """
    if A.ndim == 1:
        out = numpy.empty(A.shape[0] - 1, dtype=A.dtype)
        out[:n] = A[:n]
        out[n:] = A[n + 1:]
        return out

    if A.ndim == 2 and A.shape[0] == A.shape[1]:
        m = A.shape[0] - 1
        out = numpy.empty((m, m), dtype=A.dtype)
        out[:n, :n] = A[:n, :n]
        out[:n, n:] = A[:n, n + 1:]
        out[n:, :n] = A[n + 1:, :n]
        out[n:, n:] = A[n + 1:, n + 1:]
        return out

    ## Any other rank or a non-square 2-D array: fall through to the general path.
    return None


def remove_row_col(matrices, n, toolkit):
    result = []
    for A in matrices:
        ## The fast path is guarded on the array being a real numpy ndarray.  A
        ## symbolic toolkit hands over sympy matrices, where slicing does not mean
        ## the same thing, so those keep the generic `toolkit.delete` path -- which
        ## is also what any future array type gets until someone measures it.
        reduced = None
        if isinstance(A, numpy.ndarray) and A.dtype != object:
            reduced = _reduce_ndarray(A, n)
        if reduced is None:
            reduced = A
            for axis in range(len(A.shape)):
                reduced = toolkit.delete(reduced, [n], axis=axis)
        result.append(reduced)
    return tuple(result)

class Analysis(sim.Analysis):
    parameters = [Parameter(name='analysis', desc='Analysis name', 
                            default=None),
                  Parameter(name='epar', desc='Environment parameters',
                            default=defaultepar),
                  Parameter(name='nrsolver',
                            desc='Newton-Raphson solver strategy (a NonLinearSolver '
                                 'instance, e.g. StandardNewton(), DampedNewton(), '
                                 'GminSteppingNewton(...)); default StandardNewton()',
                            default=None),
                  Parameter(name='linearsolver',
                            desc='Linear solver strategy (a LinearSolver instance, '
                                 'e.g. DenseSolver(), SuperLUSolver(), AutoSolver()); '
                                 'default DenseSolver(), which is the historical '
                                 'numpy.linalg.solve path',
                            default=None),
                  Parameter(name='scaler',
                            desc='Jacobian scaling strategy (a Scaler instance, e.g. '
                                 'NoneScaler(), RowMaxScaler(), RowL2Scaler(), '
                                 'SinkhornKnoppScaler()); default NoneScaler()',
                            default=None)]

    def __init__(self, cir, toolkit=None, **kvargs):
        
        self.parameters = super().parameters + self.parameters
        super().__init__(cir, **kvargs)

        if toolkit is None:
            if cir.toolkit is None:
                toolkit = numeric
            else:
                toolkit = cir.toolkit

        self.toolkit = toolkit

        epar = self.par.epar
        if epar is defaultepar:
            epar = epar.copy()
            self.par.epar = epar
            
        try:
            bypass = self.par.bypass
        except (AttributeError, KeyError):
            bypass = False

        try:
            bypasstol = self.par.bypasstol
        except (AttributeError, KeyError):
            bypasstol = None
            
        if not bypass:
            bypasstol = -1.0
        elif bypasstol is None:
            # Dynamically derive bypasstol from reltol if not explicitly set
            try:
                reltol = self.par.reltol
            except (AttributeError, KeyError):
                reltol = 1e-4
            bypasstol = reltol * 1e-8
            
        epar.bypasstol = bypasstol

        if hasattr(toolkit, 'setup_analysis'):
            toolkit.setup_analysis(epar)

        self.cir = cir
        self.result = None
        self.epar = epar

    def _get_nrsolver(self):
        from pycircuit.circuit.nrsolver import NonLinearSolver, StandardNewton
        solver = getattr(self.par, 'nrsolver', None)
        if solver is None:
            return StandardNewton()
        if not isinstance(solver, NonLinearSolver):
            raise TypeError(
                "nrsolver must be a NonLinearSolver instance (e.g. StandardNewton(), "
                "DampedNewton(), GminSteppingNewton(...)), not %r" % (solver,))
        return solver

    def _get_linearsolver(self):
        from pycircuit.circuit.linearsolver import LinearSolver, DenseSolver
        solver = getattr(self.par, 'linearsolver', None)
        if solver is None:
            ## DenseSolver, not AutoSolver: the default must not change any
            ## existing result.  `numpy.linalg.solve` and SuperLU round
            ## differently, so selecting sparse automatically would move the last
            ## bits of every circuit large and sparse enough to qualify.  Opting
            ## in is the caller's decision -- see stage 7b.
            return DenseSolver()
        if not isinstance(solver, LinearSolver):
            raise TypeError(
                "linearsolver must be a LinearSolver instance (e.g. DenseSolver(), "
                "SuperLUSolver(), AutoSolver()), not %r" % (solver,))
        return solver

    def _get_scaler(self):
        from pycircuit.circuit.scaler import Scaler, NoneScaler
        scaler = getattr(self.par, 'scaler', None)
        if scaler is None:
            return NoneScaler()
        if not isinstance(scaler, Scaler):
            raise TypeError(
                "scaler must be a Scaler instance (e.g. NoneScaler(), RowMaxScaler(), "
                "RowL2Scaler(), SinkhornKnoppScaler()), not %r" % (scaler,))
        return scaler

def fsolve(f, x0, args=(), full_output=False, maxiter=200,
           xtol=1e-6, reltol=1e-4, abstol=1e-12, toolkit='Numeric', limiter=None):
    """Solve a multidimensional non-linear equation with Newton-Raphson's method

    In each iteration the linear system

    M{J(x_n)(x_{n+1}-x_n) + F(xn) = 0

    is solved and a new value for x is obtained x_{n+1}
    
    """
    
    converged = False
    ier = 2
    for i in range(maxiter):
        F, J = f(x0, *args) # TODO: Make sure J is never 0, e.g. by gmin (stepping)
        xdiff = toolkit.linearsolver(J, -F)# TODO: Limit xdiff to improve convergence

        x = x0 + xdiff
        
        if limiter is not None:
            x = limiter(x, x0)

        # KCL Scale: Upper bound of absolute branch currents/voltages
        I_scale = toolkit.dot(abs(J), abs(x)) + abs(F)

        conv_x = toolkit.alltrue(abs(xdiff) < reltol * toolkit.maximum(abs(x), abs(x0)) + xtol)
        conv_f = toolkit.alltrue(abs(F) < reltol * I_scale + abstol)

        if conv_x and conv_f:
            ier = 1
            mesg = "Success"
            break
            
        x0 = x

    if ier == 2:
        mesg = "No convergence. xerror = "+str(xdiff)
    
    infodict = {}
    if full_output:
        return x, infodict, ier, mesg
    else:
        return x


if __name__ == "__main__":
    import doctest
    doctest.testmod()

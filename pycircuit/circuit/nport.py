# -*- coding: latin-1 -*-
# Copyright (c) 2008 Pycircuit Development Team
# See LICENSE for details.

from __future__ import division

import numpy as np
from copy import copy
from . import constants as constants


def _inv(M):
    """Matrix inverse of ``M``, which may hold sympy expressions.

    ``np.linalg.inv`` requires a float/complex dtype; an object-dtype array
    (a symbolic toolkit's S-parameters) needs ``sympy.Matrix.inv()`` instead.
    """
    M = np.asarray(M)
    if M.dtype == object:
        import sympy
        return np.array(sympy.Matrix(M.tolist()).inv().tolist(), dtype=object)
    return np.linalg.inv(M)


def _is_swept(M):
    """True when ``M`` holds Waveform entries -- a frequency-swept result.

    ⚠ A swept `solve_s(freqs=array)` stores one Waveform per matrix entry, in
    an object array.  `np.linalg.inv` cannot take that and `_inv` hands it to
    sympy, which cannot either -- so every conversion has to go frequency by
    frequency (2026-09-15; `.Z` and `.Y` both raised on a sweep before).
    """
    M = np.asarray(M)
    if M.dtype != object or M.size == 0:
        return False
    return all(hasattr(x, 'get_x') and hasattr(x, 'get_y') for x in M.ravel())


def _per_frequency(M, convert):
    """Apply ``convert`` (a single-frequency matrix map) at every frequency of
    a Waveform-valued matrix, and return Waveforms on the SAME axis, with the
    source's labels and units."""
    from pycircuit.post.waveform import Waveform
    M = np.asarray(M)
    probe = M.flat[0]
    x = probe.get_x()
    nf = len(np.asarray(probe.get_y()))
    rows, cols = M.shape
    cube = np.empty((nf, rows, cols), dtype=complex)
    for a in range(rows):
        for b in range(cols):
            cube[:, a, b] = np.asarray(M[a, b].get_y(), dtype=complex)
    out = np.array([np.asarray(convert(cube[k]), dtype=complex)
                    for k in range(nf)])
    res = np.empty(out.shape[1:], dtype=object)
    for a in range(res.shape[0]):
        for b in range(res.shape[1]):
            res[a, b] = Waveform(x=x, y=out[:, a, b], xlabels=probe.xlabels,
                                 xunits=probe.xunits)
    return res


class NPort(object):
    """Class that represents an n-port with optional noise parameters

    Attributes
    ----------
    n -- number of ports
    passive -- True if n-port is passive
    noise -- True if n-port has noise parameters 
    S -- S-parameter matrix
    Y -- Y-parameter matrix
    Z -- Z-parameter matrix
    A -- ABCD-parameter matrix
    CS -- Noise wave correlation matrix
    CY -- Y-parameter noise correlation matrix
    CZ -- Z-parameter noise correlation matrix
    CA -- ABCD-parameter noise correlation matrix

    """
    passive = False

    def __init__(self, passive = False):
        self.passive = passive
    
    def __mul__(self, a):
        """Cascade of two n-ports"""
        selfA = np.asarray(self.A)
        aA = np.asarray(a.A)
        anport = NPortA(selfA @ aA,
                        selfA @ np.asarray(a.CA) @ selfA.conj().T + self.CA )
        return self.__class__(anport)

    def __floordiv__(self, a):
        """Parallel of two n-ports

        >>> import sympy as S
        >>> a,b,c,d = S.symbols('abcd')
        >>> A = NPortA(np.array([[a,b], [c,d]]))
        >>> A = S.Matrix((A // A).A)
        >>> A.simplify()
        >>> A
        [    a, 0.5*b]
        [2.0*c,     d]

        """
        ynport = NPortY(self.Y + a.Y, self.CY + a.CY)
        return self.__class__(ynport)

    def series(self, a):
        """Series connection with another n-port

        >>> import sympy as S
        >>> a,b,c,d = S.symbols('abcd')
        >>> A = NPortA(np.array([[a,b], [c,d]]))
        >>> A = S.Matrix(A.series(A).A)
        >>> A.expand()
        [    a, 2*b]
        [0.5*c,   d]
        
        """
        znport = NPortZ(self.Z + a.Z, self.CZ + a.CZ)
        
        return self.__class__(znport)

    def noisy_passive_nport(self, T=290):
        """Returns an n-port with noise parameters set if passive"""

        if not self.passive:
            raise ValueError("Cannot calculate noise-correlation matrix of "
                             "non-passive n-port")
            
        ynport = NPortY(self)
        
        ynport.CY = 4 * constants.kboltzmann * T * np.real(ynport.Y)

        return ynport        

class NPortY(NPort):
    """Two-port class where the internal representation is the Y-parameters"""
    
    def __init__(self, Y, CY = None, passive = False):
        self.passive = passive
        if isinstance(Y, NPort):
            self.Y = Y.Y
            self.CY = Y.CY
        else:
            self.Y = np.array(Y)
        
            if CY is None:
                self.CY = np.zeros(np.shape(self.Y))
            else:
                self.CY = np.array(CY)

        self.n = np.size(self.Y,0)
        
    @property
    def A(self):
        """Return chain parameters (ABCD)"""
        if self.n != 2:
            raise ValueError('N-port must be a 2-port')

        Y = self.Y
        d = Y[0,0] * Y[1,1] - Y[0,1] * Y[1,0]
        return np.array([[-Y[1,1] / Y[1,0], -1.0 / Y[1,0]],
                          [-d / Y[1,0], -Y[0,0] / Y[1,0]]])

    @property
    def Z(self):
        """Return Z-parameter matrix"""
        return np.linalg.inv(self.Y)

    @property
    def S(self, z0 = 50.0):
        """Return scattering parameters"""
        Y = self.Y
        E = np.eye(self.n, self.n)
        Zref = z0 * E
        Gref = 1 / np.sqrt(np.real(z0)) * E
        return Gref @ (E - Zref @ Y) @ np.linalg.inv(E + Zref @ Y) @ \
            np.linalg.inv(Gref)

    @property
    def CZ(self):
        Z = np.asarray(self.Z)
        return np.asarray(Z @ self.CY @ Z.conj().T)

    @property
    def CS(self, z0 = 50.):
        S = np.asarray(self.S)
        E = np.eye(self.n, self.n)
        return np.asarray((E + S) @ (np.asarray(self.CY) * z0) @
                          (E + S).conj().T / 4)

    @property
    def CA(self):
        T = np.array(self.A, copy=True)
        T[0,0] = 0
        T[1,0] = 1

        return np.asarray(T @ np.asarray(self.CY) @ T.conj().T)

class NPortZ(NPort):
    """Two-port class where the internal representation is the Z-parameters"""
    
    def __init__(self, Z, CZ = None, passive=False):
        self.passive = passive
        if isinstance(Z, NPort):
            self.Z = Z.Z
            self.CZ = np.array(Z.CZ)
        else:
            self.Z = np.array(Z)
        
            if CZ is None:
                self.CZ = np.zeros(np.shape(self.Y))
            else:
                self.CZ = np.array(CZ)

        self.n = np.size(self.Z,0)

    @property
    def A(self):
        """Return chain parameters (ABCD)"""
        if self.n != 2:
            raise ValueError('N-port must be a 2-port')

        Z = self.Z
        d = Z[0,0] * Z[1,1] - Z[0,1] * Z[1,0]
        return np.array([[Z[0,0] / Z[1,0], d / Z[1,0]],
                          [1.0 / Z[1,0], Z[1,1] / Z[1,0]]])
    
    @property
    def Y(self):
        """Return Z-parameter matrix"""
        return np.linalg.inv(self.Z)

    @property
    def S(self, z0 = 50.0):
        """Return scattering parameters"""
        Z = self.Z
        E = np.eye(self.n, self.n)
        Zref = z0 * E
        Gref = 1 / np.sqrt(np.real(z0)) * E
        return Gref @ (Z - Zref) @ np.linalg.inv(Z + Zref) @ np.linalg.inv(Gref)

    @property
    def CY(self):
        Y = np.asarray(self.Y)
        return np.asarray(Y @ self.CZ @ Y.conj().T)

    @property
    def CS(self, z0 = 50.):
        S = np.asarray(self.S)
        E = np.eye(self.n, self.n)
        T = (E - S) / (2 * np.sqrt(z0))
        return np.asarray(T @ np.asarray(self.CZ) @ T.conj().T)

    @property
    def CA(self):
       T = np.array([[1, -self.A[0,0]], [0, -self.A[1,0]]])
       return np.asarray(T @ np.asarray(self.CZ) @ T.conj().T)

class NPortA(NPort):
    """Two-port class where the internal representation is the ABCD-parameters"""

    def __init__(self, A, CA = None, passive=False):
        self.passive = passive

        if isinstance(A, NPort):
            self.A = A.A
            self.CA = A.CA
        else:
            self.A = np.array(A)

            if CA is None:
                self.CA = np.zeros(np.shape(self.Y))
            else:
                self.CA = np.array(CA)

        if np.shape(self.A) != (2,2):
            raise ValueError('Can only create ABCD-two ports')
        

        self.n = 2

    @property
    def Z(self):
        """Return Z-parameter matrix"""
        A = self.A
        d = A[0,0] * A[1,1] - A[0,1] * A[1,0]
        return np.array([[A[0,0] / A[1,0], d / A[1,0]],
                        [1.0 / A[1,0], A[1,1] / A[1,0]]])
    

    @property
    def Y(self):
        """Return Y-parameter matrix"""
        A = self.A
        d = A[0,0] * A[1,1] - A[0,1]*A[1,0]

        return np.array([[A[1,1] / A[0,1], -d / A[0,1]],
                        [-1.0 / A[0,1], A[0,0] / A[0,1]]])
    
    @property
    def S(self, z0 = 50.0):
        """Return scattering parameters
        
        >>> abcd = np.array([[  5.90000000e-01,   8.05000000e+01], \
                              [  4.20000000e-03,   1.59000000e+00]])
        >>> P = NPortA(abcd)
        >>> P.S
        array([[ 0.1,  0.3],
               [ 0.5,  0.6]])

        >>> 
        """
        a,b,c,d = self.A[0,0], self.A[0,1], self.A[1,0], self.A[1,1]

        A = np.array([[a + b / z0 - c * z0 - d, 2 * (a * d - b * c)],
                       [2,                       -a+b/z0-c*z0+d]])
        return 1/(a + b / z0 + c * z0 + d) * A

    @property
    def CY(self):
        Y = self.Y
        T = np.array([[-Y[0,0], 1], [-Y[1,0], 0]])

        return np.asarray(T @ np.asarray(self.CA) @ T.conj().T)

    @property
    def CZ(self):
        Z = np.asarray(self.Z)
        T = np.array([[1, -Z[0,0]], [0, -Z[1,0]]])
        return np.asarray(T @ np.asarray(self.CA) @ T.conj().T)

    @property
    def CS(self, z0=50.):
        return NPortY(self).CS

    def __str__(self):
        return self.__class__.__name__ + '(' + repr(self.A) + ')'


class NPortS(NPort):
    """Two-port class where the internal representation is the S-parameters"""

    def __init__(self, S, CS = None, z0 = 50, passive=False, toolkit=None):
        self.passive = passive

        self.z0 = z0

        if toolkit is None:
            from .toolkit import numeric
            toolkit = numeric
        self.toolkit = toolkit
        
        if isinstance(S, NPort):
            self.S = S.S
            self.CS = S.CS
        else:
            self.S = np.array(S)
        
            if CS is None:
                self.CS = np.zeros(np.shape(self.S))
            else:
                self.CS = np.array(CS)

        self.n = np.size(self.S,0)

    @property
    def A(self):
        """Return chain parameters (ABCD)
        
        >>> S = np.array([[0.1,0.3],[0.5,0.6]])
        >>> NPortS(S).A
        array([[0.59, 80.5],
               [0.0042, 1.59]], dtype=object)

        """
        s = self.S
        z0 = self.z0
        
        a = ((1 + s[0,0]) * (1 - s[1,1]) + s[0,1]*s[1,0]) / (2 * s[1,0])
        b = z0 * ((1 + s[0,0]) * (1 + s[1,1]) - s[0,1]*s[1,0]) / (2 * s[1,0])
        c = 1 / z0 * ((1 - s[0,0]) * (1 - s[1,1]) - s[0,1]*s[1,0]) / (2 * s[1,0])
        d = ((1 - s[0,0]) * (1 + s[1,1]) + s[0,1]*s[1,0]) / (2 * s[1,0])
        
        return np.array([[a,b],[c,d]], object)

    @property
    def Z(self):
        """Return Z-parameter matrix

        ⚠⚠ THIS CAST S TO FLOAT until 2026-09-15 (peer report, reproduced):
        `np.asarray(self.S).astype(float)` dropped Im S with only a
        ComplexWarning -- 1.027 off on a synthetic complex S, and 5103x off at
        1e7 Hz on an ordinary R-C two-port from `TwoPortAnalysis`, where
        `E - S` is nearly singular and an Im S of 3.8e-4 dominates the
        inverse.  `Y` below never had the cast.  Now written the same way as
        `Y` (S's own dtype, `_inv`, the toolkit's `sqrt`), and a swept S is
        converted frequency by frequency.
        """
        S = np.asarray(self.S)
        if _is_swept(S):
            return _per_frequency(
                S, lambda s: NPortS(s, z0=self.z0, toolkit=self.toolkit).Z)
        E = np.eye(self.n, self.n)
        zref_scalar = self.z0
        gref_scalar = 1 / self.toolkit.sqrt(self.toolkit.real(self.z0))
        Zref = zref_scalar * E
        Gref = gref_scalar * E
        Gref_inv = (1 / gref_scalar) * E
        return np.asarray(Gref_inv @ _inv(E - S) @ (S + E) @ Zref @ Gref)

    @property
    def Y(self):
        """Return Y-parameter matrix"""
        S = np.asarray(self.S)
        ## a swept S holds Waveforms, which `_inv` would hand to sympy -- see
        ## `_is_swept`; this raised on every sweep before 2026-09-15
        if _is_swept(S):
            return _per_frequency(
                S, lambda s: NPortS(s, z0=self.z0, toolkit=self.toolkit).Y)
        E = np.eye(self.n, self.n)
        ## Gref, Zref are scalar multiples of the identity -- their inverse
        ## is just the reciprocal scalar, needing no matrix inversion at
        ## all (np.linalg.inv can't handle the object-dtype arrays a
        ## symbolic toolkit's z0 produces here anyway).
        zref_scalar = self.z0
        gref_scalar = 1 / self.toolkit.sqrt(self.toolkit.real(self.z0))
        Zref_inv = (1 / zref_scalar) * E
        Gref_inv = (1 / gref_scalar) * E
        Gref = gref_scalar * E
        return np.asarray(Gref_inv @ Zref_inv @
                          _inv(S + E) @ (E - S) @ Gref)

    @property
    def CY(self):
        Y = np.asarray(self.Y)
        y0 = 1. / self.z0
        E = np.eye(self.n, self.n)
        T = (y0 * E + Y) / np.sqrt(y0)
        return np.asarray(T @ np.asarray(self.CS) @ T.conj().T)

    @property
    def CZ(self):
        Z = self.Z
        E = np.eye(self.n, self.n)
        T = (self.z0 * E + Z) / np.sqrt(self.z0)
        return np.asarray(T @ np.asarray(self.CS) @ T.conj().T)

    @property
    def CA(self):
#        return NPortZ(self).CA
        z0 = self.z0
        A = np.asarray(self.A)
        T = np.array([[np.sqrt(z0), -(A[0,1]+A[0,0]*z0)/np.sqrt(z0)],
                        [-1/np.sqrt(z0), -(A[1,1]+A[1,0]*z0)/np.sqrt(z0)]])
        return np.asarray(T @ np.asarray(self.CS) @ T.conj().T)


if __name__ == "__main__":
    import doctest
    doctest.testmod()

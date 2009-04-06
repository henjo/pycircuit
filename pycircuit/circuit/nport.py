# -*- coding: latin-1 -*-
# Copyright (c) 2008 Pycircuit Development Team
# See LICENSE for details.

import numpy as np
from copy import copy
import constants

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
        selfA = np.mat(self.A)
        aA = np.mat(a.A)
        anport = NPortA(selfA * aA,
                        selfA * np.mat(a.CA) * selfA.H + self.CA )
        return self.__class__(anport)

    def __floordiv__(self, a):
        """Parallel of two n-ports

        >>> import sympy as S
        >>> a,b,c,d = S.symbols('abcd')
        >>> A = NPortA(np.array([[a,b], [c,d]]))
        >>> (A // A).A
        array([[a, 0.5*b],
           [-2c, d]], dtype=object)

        """
        ynport = NPortY(self.Y + a.Y, self.CY + a.CY)
        return self.__class__(ynport)

    def series(self, a):
        """Series connection with another n-port

        >>> import sympy as S
        >>> a,b,c,d = S.symbols('abcd')
        >>> A = NPortA(np.array([[a,b], [c,d]]))
        >>> A.series(A).A
        [[ a 2*b]
         [ c/2 d]]
        
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
        
            if CY == None:
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
        E = np.mat(np.eye(self.n, self.n))
        Zref = z0 * E
        Gref = 1 / np.sqrt(np.real(z0)) * E
        return Gref * (E - Zref * Y) * (E + Zref * Y)**-1 * Gref**-1

    @property
    def CZ(self):
        Z = np.mat(self.Z)
        return np.array(Z * self.CY * Z.H)

    @property
    def CS(self, z0 = 50.):
        S = np.mat(self.S)
        E = np.mat(np.eye(self.n, self.n))
        return np.array((E + S) * (np.mat(self.CY)*z0) * (E + S).H / 4)

    @property
    def CA(self):
        T = np.mat(self.A)
        T[0,0] = 0
        T[1,0] = 1

        return np.array(T * np.mat(self.CY) * T.H)

class NPortZ(NPort):
    """Two-port class where the internal representation is the Z-parameters"""
    
    def __init__(self, Z, CZ = None, passive=False):
        self.passive = passive
        if isinstance(Z, NPort):
            self.Z = Z.Z
            self.CZ = np.array(Z.CZ)
        else:
            self.Z = np.array(Z)
        
            if CZ == None:
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
        E = np.mat(np.eye(self.n, self.n))
        Zref = z0 * E
        Gref = 1 / np.sqrt(np.real(z0)) * E
        return Gref * (Z - Zref) * (Z + Zref)**-1 * Gref**-1

    @property
    def CY(self):
        Y = np.mat(self.Y)
        return Y * self.CZ * Y.H

    @property
    def CS(self, z0 = 50.):
        S = np.mat(self.S)
        E = np.mat(np.eye(self.n, self.n))
        T = (E - S) / (2 * np.sqrt(z0))
        return  T * np.mat(self.CZ) * T.H

    @property
    def CA(self):
       T = np.matrix([[1, -self.A[0,0]], [0, -self.A[1,0]]])
       return np.array(T * np.mat(self.CZ) * T.H)

class NPortA(NPort):
    """Two-port class where the internal representation is the ABCD-parameters"""

    def __init__(self, A, CA = None, passive=False):
        self.passive = passive

        if isinstance(A, NPort):
            self.A = A.A
            self.CA = A.CA
        else:
            self.A = np.array(A)

            if CA == None:
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
        T = np.matrix([[-Y[0,0], 1], [-Y[1,0], 0]])

        return np.array(T * np.mat(self.CA) * T.H)

    @property
    def CZ(self):
        Z = np.mat(self.Z)
        T = np.matrix([[1, -Z[0,0]], [0, -Z[1,0]]])
        return np.array(T * np.mat(self.CA) * T.H)

    @property
    def CS(self, z0=50.):
        return NPortY(self).CS

    def __str__(self):
        return self.__class__.__name__ + '(' + repr(self.A) + ')'


class NPortS(NPort):
    """Two-port class where the internal representation is the S-parameters"""
    
    def __init__(self, S, CS = None, z0 = 50, passive=False):
        self.passive = passive

        self.z0 = z0
        
        if isinstance(S, NPort):
            self.S = S.S
            self.CS = S.CS
        else:
            self.S = np.array(S)
        
            if CS == None:
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
               [0.0, 1.59]], dtype=object)

        """
        s = self.S
        z0 = self.z0
        
        a = ((1 + s[0,0]) * (1 - s[1,1]) + s[0,1]*s[1,0]) / (2 * s[1,0])
        b = z0 * ((1 + s[0,0]) * (1 + s[1,1]) - s[0,1]*s[1,0]) / (2 * s[1,0])
        c = 1. / z0 * ((1 - s[0,0]) * (1 - s[1,1]) - s[0,1]*s[1,0]) / (2 * s[1,0])
        d = ((1 - s[0,0]) * (1 + s[1,1]) + s[0,1]*s[1,0]) / (2 * s[1,0])
        
        return np.array([[a,b],[c,d]], object)

    @property
    def Z(self):
        """Return Z-parameter matrix"""
        S = np.mat(self.S).astype(float)
        E = np.mat(np.eye(self.n, self.n))
        Zref = self.z0 * E
        Gref = 1 / np.sqrt(np.real(self.z0)) * E
        return np.array(Gref.I * (E - S).I * (S + E) * Zref * Gref)

    @property
    def Y(self):
        """Return Z-parameter matrix"""
        S = np.mat(self.S)
        E = np.mat(np.eye(self.n, self.n))
        Zref = self.z0 * E
        Gref = 1 / np.sqrt(np.real(self.z0)) * E
        return np.array(Gref**-1 * Zref**-1 * (S + E)**-1 * (E - S) * Gref)

    @property
    def CY(self):
        Y = np.mat(self.Y)
        y0 = 1. / self.z0
        E = np.mat(np.eye(self.n, self.n))
        T = (y0 * E + Y) / np.sqrt(y0)
        return np.array(T * np.mat(self.CS) * T.H)

    @property
    def CZ(self):
        Z = self.Z
        E = np.mat(np.eye(self.n, self.n))
        T = (self.z0 * E + Z) / np.sqrt(self.z0)
        return np.array(T * np.mat(self.CS) * T.H)

    @property
    def CA(self):
        return NPortZ(self).CA
        z0 = self.z0
        A = np.mat(self.A)
        T = np.matrix([[np.sqrt(z0), -(A[0,1]+A[0,0]*z0)/np.sqrt(z0)],
                        [-1/np.sqrt(z0), -(A[1,1]+A[1,0]*z0)/np.sqrt(z0)]])
        return np.array(T * np.mat(self.CS) * T.H)


if __name__ == "__main__":
    import doctest
    doctest.testmod()

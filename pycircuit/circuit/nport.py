# Copyright (c) 2008 Pycircuit Development Team
# See LICENSE for details.

import numpy as npy
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

    def __mul__(self, a):
        """Cascade of two n-ports"""
        anport = NPortA(dot(self.A, a.A), dot(dot(self.A, a.CA), npy.conjugate(self.A).T) + self.CA )
        return self.__class__(anport)

    def __floordiv__(self, a):
        """Parallel of two n-ports

        >>> import sympy as S
        >>> a,b,c,d = S.symbols('abcd')
        >>> A = NPortA(npy.array([[a,b], [c,d]]))
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
        >>> A = NPortA(npy.array([[a,b], [c,d]]))
        >>> A.series(A).A
        [[ a 2*b]
         [ c/2 d]]
        
        """
        znport = NPortZ(self.Z + a.Z, self.CZ + a.CZ)
        
        return self.__class__(znport)

    def noisy_nport(self, T=290):
        """Returns an n-port with noise parameters set if passive"""

        if not self.passive:
            raise ValueError("Cannot calculate noise-correlation matrix of non-passive n-port")
            
        ynport = NPortY(self)
        
        ynport.CY = 2 * constants.kboltzmann * T * npy.real(ynport.Y)

        return ynport        

class NPortY(NPort):
    """Two-port class where the internal representation is the Y-parameters"""
    
    def __init__(self, Y, CY = None):
        if isinstance(Y, NPort):
            self.Y = Y.Y
            self.CY = Y.CY
        else:
            self.Y = Y
        
            if CY == None:
                self.CY = npy.zeros(npy.shape(self.Y))
            else:
                self.CY = CY

        self.n = npy.size(self.Y,0)
        
    @property
    def A(self):
        """Return chain parameters (ABCD)"""
        if self.n != 2:
            raise ValueError('N-port must be a 2-port')

        Y = self.Y
        d = Y[0,0] * Y[1,1] - Y[0,1] * Y[1,0]
        return npy.array([[-Y[1,1] / Y[1,0], -1.0 / Y[1,0]],
                          [-d / Y[1,0], -Y[0,0] / Y[1,0]]])

    @property
    def Z(self):
        """Return Z-parameter matrix"""
        return npy.linalg.inv(self.Y)

    @property
    def S(self, z0 = 50.0):
        """Return scattering parameters"""
        Y = self.Y
        E = npy.mat(npy.eye(self.n, self.n))
        Zref = z0 * E
        Gref = 1 / npy.sqrt(npy.real(z0)) * E
        return Gref * (E - Zref * Y) * (E + Zref * Y)**-1 * Gref**-1

    @property
    def CZ(self):
        Z = npy.mat(self.Z)
        return npy.array(Z * self.CY * Z.H)

    @property
    def CS(self, z0 = 50.):
        S = npy.mat(self.S)
        E = npy.mat(npy.eye(self.n, self.n))
        return npy.array((E + S) * (npy.mat(self.CY)*z0) * (E + S).H / 4)

    @property
    def CA(self):
        T = npy.mat(self.A)
        T[0,0] = 0
        T[1,0] = 1

        return npy.array(T * npy.mat(self.CY) * T.H)

class NPortZ(NPort):
    """Two-port class where the internal representation is the Z-parameters"""
    
    def __init__(self, Z, CZ = None):
        if isinstance(Z, NPort):
            self.Z = Z.Z
            self.CZ = Z.CZ
        else:
            self.Z = Z
        
            if CZ == None:
                self.CZ = npy.zeros(npy.shape(self.Y))
            else:
                self.CZ = CZ

        self.n = npy.size(self.Z,0)

    @property
    def A(self):
        """Return chain parameters (ABCD)"""
        if self.n != 2:
            raise ValueError('N-port must be a 2-port')

        Z = self.Z
        d = Z[0,0] * Z[1,1] - Z[0,1] * Z[1,0]
        return npy.array([[Z[0,0] / Z[1,0], d / Z[1,0]],
                          [1.0 / Z[1,0], Z[1,1] / Z[1,0]]])
    
    @property
    def Y(self):
        """Return Z-parameter matrix"""
        return npy.linalg.inv(self.Z)

    @property
    def S(self, z0 = 50.0):
        """Return scattering parameters"""
        Z = self.Z
        E = npy.mat(npy.eye(self.n, self.n))
        Zref = z0 * E
        Gref = 1 / npy.sqrt(npy.real(z0)) * E
        return Gref * (Z - Zref) * (Z + Zref)**-1 * Gref**-1

    @property
    def CY(self):
        Y = npy.mat(self.Y)
        return Y * self.CZ * Y.H

    @property
    def CS(self, z0 = 50.):
        S = npy.mat(self.S)
        E = npy.mat(npy.eye(self.n, self.n))
        T = (E - S) / (2 * npy.sqrt(z0))
        return  T * npy.mat(self.CZ) * T.H

    @property
    def CA(self):
       T = npy.matrix([[1, -self.A[0,0]], [0, -self.A[1,0]]])
       return T * npy.mat(self.CZ) * T.H

class NPortA(NPort):
    """Two-port class where the internal representation is the ABCD-parameters"""

    def __init__(self, A, CA = None):
        if isinstance(A, NPort):
            self.A = A.A
            self.CA = A.CA
        else:
            self.A = A

            if CA == None:
                self.CA = npy.zeros(npy.shape(self.Y))
            else:
                self.CA = CA

        if npy.shape(self.A) != (2,2):
            raise ValueError('Can only create ABCD-two ports')
        

        self.n = 2

    @property
    def Z(self):
        """Return Z-parameter matrix"""
        A = self.A
        d = A[0,0] * A[1,1] - A[0,1] * A[1,0]

        return npy.array([[A[0,0] / A[1,0], d / A[1,0]],
                        [1.0 / A[1,0], A[1,1] / A[1,0]]])
    

    @property
    def Y(self):
        """Return Y-parameter matrix"""
        A = self.A
        d = A[0,0] * A[1,1] - A[0,1]*A[1,0]

        return npy.array([[A[1,1] / A[0,1], -d / A[0,1]],
                        [-1.0 / A[0,1], A[0,0] / A[0,1]]])
    
    @property
    def S(self, z0 = 50.0):
        """Return scattering parameters
        
        >>> abcd = npy.array([[  5.90000000e-01,   8.05000000e+01], \
                              [  4.20000000e-03,   1.59000000e+00]])
        >>> P = NPortA(abcd)
        >>> P.S
        array([[ 0.1,  0.3],
               [ 0.5,  0.6]])

        >>> 
        """
        a,b,c,d = self.A[0,0], self.A[0,1], self.A[1,0], self.A[1,1]

        A = npy.array([[a + b / z0 - c * z0 - d, 2 * (a * d - b * c)],
                       [2,                       -a+b/z0-c*z0+d]])
        return 1/(a + b / z0 + c * z0 + d) * A

    @property
    def CY(self):
        Y = self.Y
        T = npy.matrix([[-Y[0,0], 1], [-Y[1,0], 0]])

        return npy.array(T * npy.mat(self.CA) * T.H)

    @property
    def CZ(self):
        Z = self.Z
        T = npy.matrix([[1, -Z[0,0]], [0, -Z[1,0]]])

        return npy.array(T * npy.mat(self.CA) * T.H)

    @property
    def CS(self, z0=50.):
        return NPortY(self).CS

    def __str__(self):
        return self.__class__.__name__ + '(' + repr(self.A) + ')'


class NPortS(NPort):
    """Two-port class where the internal representation is the S-parameters"""
    
    def __init__(self, S, CS = None, z0 = 50):
        self.z0 = z0
        
        if isinstance(S, NPort):
            self.S = S.S
            self.CS = S.CS
        else:
            self.S = S
        
            if CS == None:
                self.CS = npy.zeros(npy.shape(self.S))
            else:
                self.CS = CS

        self.n = npy.size(self.S,0)

    @property
    def A(self):
        """Return chain parameters (ABCD)
        
        >>> S = npy.array([[0.1,0.3],[0.5,0.6]])
        >>> NPortS(S).A
        array([[0.59, 80.5],
               [0.0, 1.59]], dtype=object)

        """
        s = self.S
        z0 = self.z0
        
        a = ((1 + s[0,0]) * (1 - s[1,1]) + s[0,1]*s[1,0]) / (2 * s[1,0])
        b = z0 * ((1 + s[0,0]) * (1 + s[1,1]) - s[0,1]*s[1,0]) / (2 * s[1,0])
        c = 1 / z0 * ((1 - s[0,0]) * (1 - s[1,1]) - s[0,1]*s[1,0]) / (2 * s[1,0])
        d = ((1 - s[0,0]) * (1 + s[1,1]) + s[0,1]*s[1,0]) / (2 * s[1,0])
        
        return npy.array([[a,b],[c,d]], object)

    @property
    def Z(self):
        """Return Z-parameter matrix"""
        S = npy.mat(self.S)
        E = npy.mat(npy.eye(self.n, self.n))
        Zref = self.z0 * E
        Gref = 1 / npy.sqrt(npy.real(self.z0)) * E
        return npy.array(Gref.I * (E - S).I * (S + E) * Zref * Gref)

    @property
    def Y(self):
        """Return Z-parameter matrix"""
        S = npy.mat(self.S)
        E = npy.mat(npy.eye(self.n, self.n))
        Zref = self.z0 * E
        Gref = 1 / npy.sqrt(npy.real(self.z0)) * E
        return npy.array(Gref**-1 * Zref**-1 * (S + E)**-1 * (E - S) * Gref)

    @property
    def CY(self):
        Y = npy.mat(self.Y)
        y0 = 1. / self.z0
        E = npy.mat(npy.eye(self.n, self.n))
        T = (y0 * E + Y) / npy.sqrt(y0)
        return T * npy.mat(self.CS) * T.H

    @property
    def CZ(self):
        Z = self.Z
        E = npy.mat(npy.eye(self.n, self.n))
        T = (self.z0 * E + Z) / npy.sqrt(self.z0)
        return T * npy.mat(self.CS) * T.H

    @property
    def CA(self):
        return NPortZ(self).CA
        z0 = self.z0
        A = npy.mat(self.A)
        T = npy.matrix([[npy.sqrt(z0), -(A[0,1]+A[0,0]*z0)/npy.sqrt(z0)],
                        [-1/npy.sqrt(z0), -(A[1,1]+A[1,0]*z0)/npy.sqrt(z0)]])
        return  T * npy.mat(self.CS) * T.H


if __name__ == "__main__":
    import doctest
    doctest.testmod()

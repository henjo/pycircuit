import numpy as np
from scipy import factorial, interpolate, linalg


def fourier_analysis(t, x, harmonics=range(9)):
    """Evaluate the fourier series of fitted piecewise polynomial"""
    t = np.array(t)

    M = 1
    N = np.size(x,0)
    T = t[-1] - t[0]
    omega = 2 * np.pi / T
    
    ## Pre-calculate upper and lower integration limits
    t_intervals = np.array((t[:-1], t[1:])).T

    def alpha(k,m,n):
        tau = t_intervals[n-1]
        i = np.arange(m+1)[:,np.newaxis]
        the_sum = np.sum((1j*k*omega)**(-m+i-1) * \
                             tau[:,np.newaxis]**i / factorial(i), axis=-2)
        integral = -factorial(m) * np.exp(-1j * k * omega * tau) * the_sum
        return integral[:,1] - integral[:,0]

    ## Time step
    h = np.diff(t)
    ## Piecewise linear fitting x(n) + (x(n+1)-x(n)) / h * t
    xdelta = np.diff(x)
    pwc = np.array([x[:-1] - t[:-1] * xdelta / h, xdelta / h])

#    tck = interpolate.splrep(t, x, per=1)  
#    if xdot == None:
#        xdot = interpolate.splev(t, tck, der=1)

    c_k = []
    freqs = []
    for k in harmonics:
        if k == 0:
            if True:
                n = np.arange(1, N)
                tmp = np.zeros(n.shape, dtype=np.complex)

                m = np.arange(M+1)[:,np.newaxis]

                tau = t_intervals[n-1][..., np.newaxis, :]

                coeffs = pwc.T[...,np.newaxis]

                ## Sum over time instants
                integral = np.sum(coeffs / (m+1) * tau ** (m + 1), axis=0)
                ## Sum over terms
                integral = np.sum(integral, axis=-2)

                c_k.append(1 / T * (integral[...,1] - integral[..., 0]))
        else:
            n = np.arange(1, N)
            tmp = np.zeros(n.shape, dtype=np.complex)
            for m in range(M+1):
                tmp += pwc[m,n-1] * alpha(k,m,n)
            c_k.append(1/T * sum(tmp))
        
        freqs.append(1/T * k)

    return np.array(c_k), np.array(freqs)


from pycircuit.circuit import *
import numpy as np
import pylab
from pycircuit.post import from_arrays

circuit.default_toolkit = numeric

Rs=50.
Vs=1e-4


def build_lna_degen_gm():
    c = SubCircuit()
    c['Rs'] = R(1, gnd, r=Rs)
    c['vs'] = VS(2, 1, vac=Vs)
    c['vccs'] = VCCS(2, 4, 3, 4, gm=0.1)
    c['Cgs'] = C(2, 4, c=1e-9)
    c['rl'] = R(3, gnd, r=200.)
    c['rd'] = R(4, gnd, r=10.)
    return c

def wave_conv(w1,w2):
    '''Convolution of two waveforms.

    Convolution by frequency of two waveforms.
    '''
    (x1vec,), y1vec = w1.to_arrays()  # frequencies, values
    (x2vec,), y2vec = w2.to_arrays()
    wdict={}
    for x1,y1 in zip(x1vec,y1vec):
        for x2,y2 in zip(x2vec,y2vec):
            if x1+x2 in wdict:
                wdict[x1+x2]+=y1*y2
            else:
                wdict[x1+x2]=y1*y2
    newx=np.sort(np.array(list(wdict.keys())))
    newy=np.array([])
    for x in newx:
        newy=np.append(newy,wdict[x])
    newwave=from_arrays(newx,newy,xlabels=['frequency'])
    return newwave

def solve_by_freq(w,c,vac=True):
    ac=AC(c)
    (xvec,), _ = w.to_arrays()
    vc=np.array([])
    vout=np.array([])
    if vac:
        c['inl']=IS(3,4,iac=0.)
        c['vs']=VS(2,1,vac=1.)
    else:
        c['vs']=VS(2,1,vac=0.)
        c['inl']=IS(3,4,iac=1.)
    res=ac.solve(xvec)
    wvc=res.v(2,4)*w
    wvout=res.v(3)*w
    return wvc,wvout

def lna_volterra(w):
    '''Convolution based 2nd and 3rd order volterra analysis.
    
    w is a Waveform object with stimuli frequencies and amplitudes
    Using (sparse) convolution to calculate nonlinear currents.
    '''
    c=build_lna_degen_gm()
    K2gm=0.01
    K3gm=0.001
    # 1st order response of controlling voltage and output
    vc1,vout_1=solve_by_freq(w,c)
    #self convolution of 1st order spectrum
    v1_c_v1=wave_conv(vc1,vc1)
    #nonlinear current per frequency
    INL2=K2gm*v1_c_v1
    # 2nd order response of controlling voltage and output
    vc2,vout_2=solve_by_freq(INL2,c,vac=False)
    #calulate 3rd order from K3gm
    INL31=K3gm*wave_conv(vc1,v1_c_v1)
    vc31,vout_31=solve_by_freq(INL31,c,vac=False)
    #calculate 3rd order from K2gm
    INL32=K2gm*wave_conv(vc2,vc1)
    vc32,vout_32=solve_by_freq(INL32,c,vac=False)
    return vout_1,vout_2,vout_31,vout_32

def test_volterra():
    f=np.array([-1200.,-1000.,1000.,1200.])
    v=np.array([2.,2.,2.,2.])
    w=from_arrays(f,v,xlabels=['frequency'])
    v1,v2,v31,v32=lna_volterra(w)
    abs(v1).stem(ax=pylab.subplot(4,1,1))
    abs(v2).stem(ax=pylab.subplot(4,1,2))
    abs(v31).stem(ax=pylab.subplot(4,1,3))
    abs(v32).stem(ax=pylab.subplot(4,1,4))
    pylab.show()
    
if __name__ == '__main__':
    test_volterra()

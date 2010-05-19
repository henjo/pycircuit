#from numpy import loadtxt
# import I of V as array with V=A[:,0], I=A[:,1], Insq=A[:,2]
from pycircuit.circuit import *
from pycircuit.circuit import func

class myVCCS(Circuit):
   """Example of VCCS using lookup table

   >>> import pylab
   >>> import np
   >>> from pycircuit.circuit.elements import *
   >>> from pycircuit.post import plotall
   >>> from pycircuit.circuit.transient import Transient
   >>> vvec=np.linspace(-2,2,100)
   >>> ivec=np.tanh(vvec)
   >>> nvec=ivec*0 # no noise
   >>> c = SubCircuit()
   >>> n1,n2 = c.add_nodes('1', '2')
   >>> c['vsin'] = VSin(n1, gnd, freq=2e3, va=1, vo=1)
   >>> c['vccs'] = myVCCS(n1, gnd, n2, gnd, ivec=ivec, vvec=vvec, nvec=nvec)
   >>> c['rl'] = R(n2, gnd, r=2.0)
   >>> tran = Transient(c)
   >>> res = tran.solve(tend=1e-3, timestep=1e-5)
   >>> plotall(res.v(n1),res.v(n2))
   >>> pylab.show()
   
   """
   terminals = ('inp', 'inn', 'outp', 'outn')
   instparams = [Parameter(name='ivec', desc='Current vector',
                           unit='A', default=1e-3),
                 Parameter(name='vvec', desc='Dependent voltage variable vector',
                           unit='V', default=0),
                 Parameter(name='nvec',
                           desc='Current noise power spectral density vector',
                           unit='A^2/Hz', default=0.0)]
   linear = False
       
   def update(self, subject):
      self.function = func.TabFunc(self.ipar.vvec, self.ipar.ivec)
      self.noiseFunction = func.TabFunc(self.ipar.vvec, self.ipar.nvec)
      n = self.n
      G = self.toolkit.zeros((n,n))
      inpindex, innindex, outpindex, outnindex = \
          (self.nodes.index(self.nodenames[name]) 
           for name in ('inp', 'inn', 'outp', 'outn'))
      G[outpindex, inpindex] += 1
      G[outpindex, innindex] += -1
      G[outnindex, inpindex] += -1
      G[outnindex, innindex] += 1
      self._G = G
      CY = self.toolkit.zeros((n,n))
      CY[outpindex, outpindex] +=1
      CY[outpindex, outnindex] +=-1
      CY[outnindex, outnindex] +=1
      CY[outnindex, outpindex] +=-1
      self._CY = CY
      self._I = self.toolkit.zeros(n)
      self._I[outpindex] +=-1
      self._I[outnindex] += 1

   def G(self, x, epar=defaultepar):
      gm=self.function.fprime(x[1]-x[0])
      return self._G*gm

   def K2G(self, x, epar=defaultepar):
      K2gm=self.function.fprime(x[1]-x[0],order=2)
      return self._G*K2gm

   def K3G(self, x, epar=defaultepar):
      K3gm=self.function.fprime(x[1]-x[0],order=3)
      return self._G*K3gm

   def i(self, x, epar=defaultepar):
      """
      
      """
      i=self.function.f(x[1]-x[0])
      return self._I*i

   def CY(self, x, w, epar=defaultepar):
      #gm=self.function.fprime(x[1]-x[0])
      #ipsd=4*kT*gm
      ipsd=self.noiseFunction.f(x[1]-x[0])
      return  self._CY*ipsd


if __name__ == "__main__":
    import doctest
    doctest.testmod()

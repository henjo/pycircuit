# Copyright (c) 2008 Pycircuit Development Team
# See LICENSE for details.

import unittest
import copy

class Parameter(object):
    def __init__(self, name, desc=None, unit=None, default=None):
        self.name=name
        self.desc=desc
        self.default=default
        self.unit=unit

    def __hash__(self):
        return self.name.__hash__()

    def copy(self):
        return copy.copy(self)

    def __repr__(self):
        return self.name

class ParameterDict(object):
    def __init__(self, *parameters, **kvargs):
        self._parameters = {}
        self._paramnames = []
        self.append(*parameters)
        self.set(**kvargs)
        
        
    def append(self, *parameters):
        for param in parameters:
            if param.name not in self._parameters:
                self._paramnames.append(param.name)
                self._parameters[param.name] = param
                self.__dict__[param.name] = param.default
                
    def set(self, **kvargs):
        for k,v in kvargs.items():
            if k not in self.__dict__:
                raise KeyError('parameter %s not in parameter dictionary'%k )
            self.__dict__[k] = v

    def get(self, param):
        """Get value by parameter object or parameter name"""
        if isinstance(param, Parameter):
            return self.__dict__[param.name]
        else:
            return self.__dict__[param]

    def copy(self, *parameters, **kvargs):
        newpd = ParameterDict()
        newpd.__dict__ = copy.copy(self.__dict__)
        newpd._parameters = copy.copy(self._parameters)
        newpd._paramnames = copy.copy(self._paramnames)
        newpd.append(*parameters)
        newpd.set(**kvargs)

        print newpd._parameters is self._parameters
        return newpd

    def __getitem__(self, key):
        return self._parameters[key]

    def __contains__(self, key):
        if isinstance(key, Parameter):
            return key.name in self._parameters and self._parameters[key.name] == key 
        return key in self._parameters
    
    def __len__(self):
        return len(self._paramnames)

    parameters = property(lambda self: [self._parameters[name] for name in self._paramnames])

class ParameterTest(unittest.TestCase):
    """Test Parameter class"""

    def testInstance(self):
        """Testing instantiating Parameter objects"""
        param1 = Parameter("gm", "Transconductance", "A/V", default=1e-3)
        param2 = Parameter(name="gm", desc="Transconductance", unit="A/V", default=1e-3)
        for param in param1,param2:
            self.assertEqual(param.name, 'gm')
            self.assertAlmostEqual(param.default, 1e-3)
            self.assertEqual(param.desc, 'Transconductance')
            self.assertEqual(param.unit, 'A/V')

    def testIncorrectInstantiation(self):
        def testIncorrect():
            param = Parameter(apa=1)

        self.failUnlessRaises(TypeError, testIncorrect)
        
class ParameterDictTest(unittest.TestCase):
    def testInstantiate(self):
        paramdict = ParameterDict()
        gm = Parameter(name="gm", desc="Transconductance", unit="A/V", default=1e-3)
        gds = Parameter(name="gds", desc="Output conductance", unit="A/V", default=1e-6)
        paramdict = ParameterDict(gm,gds,
                                  gm=2e-3
                                  )
        self.assertEqual(paramdict.parameters, [gm,gds])
        self.assertEqual(paramdict.gm, 2e-3)
        self.assertEqual(paramdict.gds, 1e-6)

    def testAppendParameter(self):
        """Test appending a parameter"""
        paramdict = ParameterDict()
        self.assertEqual(len(paramdict), 0)
        gmparam = Parameter(name="gm", desc="Transconductance", unit="A/V")
        paramdict.append(gmparam)
        self.assertEqual(len(paramdict), 1)
        self.assertEqual(paramdict['gm'], gmparam)

    def testGetItem(self):
        """Test getting parameter"""
        paramdict = ParameterDict()
        gmparam = Parameter(name="gm", desc="Transconductance", unit="A/V", default=2.0e-6)
        paramdict.append(gmparam)
        self.assertEqual(paramdict['gm'], gmparam)

    def testGetSetValue(self):
        paramdict = ParameterDict()
        gmparam = Parameter(name="gm", desc="Transconductance", unit="A/V", default=2.0e-6)
        paramdict.append(gmparam)

        self.assertEqual(paramdict.get('gm'), 2e-6)

        self.assertEqual(paramdict.get(gmparam), 2e-6)

        self.assertEqual(paramdict.gm, 2e-6)

        paramdict.gm = 3e-6
        self.assertEqual(paramdict.gm, 3e-6)

    def testSet(self):
        paramdict = ParameterDict()
        paramdict.append(Parameter(name="gm", desc="Transconductance", unit="A/V", default=2.0e-6))
        paramdict.set(gm=4e-6)
        self.assertEqual(paramdict.gm, 4e-6)
        def test():
            paramdict.set(gds=3e-6)
        self.failUnlessRaises(KeyError, test)
        
    def testContains(self):
        paramdict = ParameterDict()
        gmparam = Parameter(name="gm", desc="Transconductance", unit="A/V", default=2.0e-6)
        paramdict.append(gmparam)
        self.assertTrue('gm' in paramdict)
        self.assertTrue(gmparam in paramdict)
        self.assertFalse('gds' in paramdict)
        self.assertFalse(Parameter(name='gds') in paramdict)

    def testGetParameters(self):
        paramdict = ParameterDict()
        gmparam = Parameter(name="gm", desc="Transconductance", unit="A/V", default=2.0e-6)
        gdsparam = Parameter(name="gds")
        paramdict.append(gmparam)
        self.assertEqual(paramdict.parameters, [gmparam])
        paramdict.append(gdsparam)
        self.assertEqual(paramdict.parameters, [gmparam, gdsparam])

    def testCopy(self):
        paramdict = ParameterDict(Parameter(name="vth0", desc="Threshold voltage", unit="V", default=0.3))
        paramdict2 = paramdict.copy()
        paramdict2.vth0 = 0.4
        paramdict3 = paramdict2.copy(vth0=0.5)
        self.assertEqual(paramdict2.vth0, 0.4)
        self.assertEqual(paramdict.vth0, 0.3)
        self.assertEqual(paramdict3.vth0, 0.5)

if __name__ == "__main__":
    unittest.main()



